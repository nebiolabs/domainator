"""Persisted record-offset indexes (``.didx``) and shard resolution for
domain_search databases.

A *domain_search database* is a GenBank or FASTA file, optionally BGZF-compressed,
optionally split into numbered *shards*. This module owns two things:

1. The ``.didx`` sidecar format: a small binary file holding the per-record
   ``(offset, cds_count)`` pairs that ``utils.get_offsets`` would otherwise have to
   recompute by scanning the whole database on every run. Offsets are plain byte
   offsets for uncompressed files and BGZF virtual offsets (``block_start<<16 |
   within``) for BGZF files; the index records which, so a reader can't confuse them.

2. Shard naming/resolution: a shard inserts ``.N`` after the base name, before the
   format/compression suffixes (``mydb.gb`` -> ``mydb.0.gb`` / ``mydb.0.gb.bgz``);
   leading zeros are allowed but ignored. An index appends ``.didx`` to the *full*
   file/shard name (``mydb.0.gb.bgz.didx``).

Correctness invariant: a stale or unreadable index is NEVER trusted -- a wrong
offset would make a worker seek into the middle of a record and silently corrupt
results. Readers validate a size+mtime fingerprint and fall back to a fresh scan
(with a one-time warning) on any mismatch. The format is intentionally
dependency-light and does not import :mod:`domainator.utils` (callers pass the
already-detected ``filetype`` and ``compression``), to avoid an import cycle.
"""
import os
import re
import struct
import sys
import warnings
from array import array

# Bump the major version on any breaking layout change (mirrors data_matrix
# ._MATRIX_FILE_VERSION discipline). A reader rejects a different major version and
# recomputes; a newer minor with the same major stays readable (only additive
# changes via the reserved/flags fields are allowed within a major version).
_INDEX_FILE_VERSION = "1.0"
_VER_MAJOR, _VER_MINOR = (int(x) for x in _INDEX_FILE_VERSION.split("."))

INDEX_EXTENSION = ".didx"
_MAGIC = b"DIDX"

# 56-byte header: magic, version(major,minor), filetype flag, compression flag,
# 6 reserved/pad bytes (keeps the body 8-byte aligned), source size + mtime_ns
# fingerprint, record count, total CDS/protein count (sum of the body's second
# column, so -Z 0 can read the target count without scanning the body), and a
# reserved flags bitfield. Then record_count (offset, cds_count) pairs.
_HEADER_STRUCT = struct.Struct("<4sHHBB6xQQQQQ")
_PAIR_STRUCT = struct.Struct("<QQ")
assert _HEADER_STRUCT.size == 56

_FILETYPE_TO_FLAG = {"genbank": 0, "fasta": 1}
_FLAG_TO_FILETYPE = {0: "genbank", 1: "fasta"}

# Mirrors utils.COMPRESSION_EXTENSIONS; duplicated here to keep this module free of
# a utils import (the values are part of the on-disk naming contract and stable).
_COMPRESSION_EXTENSIONS = frozenset({"gz", "bgz", "bgzf"})

# Index paths already warned about this process, so a stale/bad index warns once
# (not once per worker task). Workers in separate processes each warn once, which
# is acceptable and low-volume.
_warned_index_paths = set()


def _compression_flag(compression):
    """Map a compression label to the on-disk flag. Plain gzip is unsupported."""
    if compression == "bgzf":
        return 1
    if compression in (None, ""):
        return 0
    raise ValueError(
        f"Cannot build a .didx index for {compression!r}-compressed input; only "
        "uncompressed or BGZF (block-gzip) databases support random-access offsets."
    )


def index_path_for(source_path):
    """The ``.didx`` sidecar path for a database file/shard (append to full name)."""
    return str(source_path) + INDEX_EXTENSION


def source_fingerprint(source_path):
    """Return ``(size, mtime_ns)`` used to detect a stale index. Single ``os.stat``."""
    st = os.stat(source_path)
    return st.st_size, st.st_mtime_ns


def _warn_once(index_path, message):
    if index_path not in _warned_index_paths:
        _warned_index_paths.add(index_path)
        warnings.warn(message, RuntimeWarning, stacklevel=3)


def write_index(source_path, offsets, num_proteins, *, filetype, compression):
    """Write a ``.didx`` index for ``source_path`` and return its path.

    ``offsets``/``num_proteins`` are parallel sequences (typically ``array('Q')``)
    of record byte/virtual offsets and per-record CDS/protein counts, exactly as
    produced by ``utils.get_offsets``. ``compression`` is ``None`` or ``"bgzf"``
    (raises ``ValueError`` for plain gzip).
    """
    comp_flag = _compression_flag(compression)
    if filetype not in _FILETYPE_TO_FLAG:
        raise ValueError(f"Cannot index filetype {filetype!r}; expected genbank or fasta.")
    if len(offsets) != len(num_proteins):
        raise ValueError("offsets and num_proteins must have equal length")

    size, mtime_ns = source_fingerprint(source_path)
    n = len(offsets)
    total_cds = 0
    body = array("Q")
    for off, cds in zip(offsets, num_proteins):
        body.append(int(off))
        body.append(int(cds))
        total_cds += int(cds)
    if sys.byteorder != "little":
        body.byteswap()

    index_path = index_path_for(source_path)
    header = _HEADER_STRUCT.pack(
        _MAGIC, _VER_MAJOR, _VER_MINOR, _FILETYPE_TO_FLAG[filetype], comp_flag,
        size, mtime_ns, n, total_cds, 0,
    )
    with open(index_path, "wb") as f:
        f.write(header)
        f.write(body.tobytes())
    return index_path


def _validate_header(index_path, source_path, *, filetype, compression):
    """Open and validate the index header against the source. Return ``(handle,
    record_count, total_cds)`` if valid+fresh, else ``None`` (warning once on
    stale/corrupt). The caller owns closing the returned handle."""
    try:
        f = open(index_path, "rb")
    except OSError:
        return None  # no index -> silent miss (normal fast-path absence)
    try:
        header = f.read(_HEADER_STRUCT.size)
        if len(header) < _HEADER_STRUCT.size:
            _warn_once(index_path, f"Ignoring truncated index {index_path}; recomputing offsets.")
            f.close()
            return None
        magic, maj, _minor, ft_flag, comp_flag, size, mtime_ns, n, total_cds, _flags = _HEADER_STRUCT.unpack(header)
        if magic != _MAGIC or maj != _VER_MAJOR:
            _warn_once(index_path, f"Ignoring unreadable/incompatible index {index_path}; recomputing offsets.")
            f.close()
            return None
        if ft_flag != _FILETYPE_TO_FLAG.get(filetype) or comp_flag != _compression_flag(compression):
            _warn_once(index_path, f"Ignoring index {index_path} built for a different file kind; recomputing offsets.")
            f.close()
            return None
        # Truncation guard: file must be exactly header + n pairs.
        expected = _HEADER_STRUCT.size + n * _PAIR_STRUCT.size
        if os.fstat(f.fileno()).st_size != expected:
            _warn_once(index_path, f"Ignoring corrupt index {index_path} (wrong size); recomputing offsets.")
            f.close()
            return None
        # Staleness: the source must be unchanged since the index was built.
        try:
            cur_size, cur_mtime = source_fingerprint(source_path)
        except OSError:
            f.close()
            return None
        if (cur_size, cur_mtime) != (size, mtime_ns):
            _warn_once(index_path, f"Ignoring stale index {index_path} (source changed since it was built); recomputing offsets.")
            f.close()
            return None
        return f, n, total_cds
    except Exception:
        # Never let a malformed index raise into the offset path; recompute.
        _warn_once(index_path, f"Ignoring unreadable index {index_path}; recomputing offsets.")
        try:
            f.close()
        except Exception:
            pass
        return None


def read_index(source_path, *, filetype, compression):
    """Return ``(offsets, num_proteins)`` as two ``array('Q')`` from the valid,
    fresh ``.didx`` for ``source_path``, or ``None`` to signal the caller to scan."""
    validated = _validate_header(index_path_for(source_path), source_path,
                                 filetype=filetype, compression=compression)
    if validated is None:
        return None
    f, n, _total_cds = validated
    try:
        body = f.read(n * _PAIR_STRUCT.size)
    finally:
        f.close()
    interleaved = array("Q")
    interleaved.frombytes(body)
    if sys.byteorder != "little":
        interleaved.byteswap()
    return array("Q", interleaved[0::2]), array("Q", interleaved[1::2])


def i_read_index(source_path, *, filetype, compression):
    """Like :func:`read_index` but stream ``(offset, cds_count)`` tuples lazily
    (memory-light for huge databases). Returns ``None`` on miss/stale/corrupt."""
    validated = _validate_header(index_path_for(source_path), source_path,
                                 filetype=filetype, compression=compression)
    if validated is None:
        return None
    f, n, _total_cds = validated

    def _gen():
        try:
            chunk_pairs = 1 << 16
            chunk_bytes = chunk_pairs * _PAIR_STRUCT.size
            remaining = n
            while remaining > 0:
                want = min(remaining, chunk_pairs) * _PAIR_STRUCT.size
                buf = f.read(want)
                if len(buf) < want:
                    return  # truncated mid-stream; stop (header size-check should prevent)
                for off, cds in _PAIR_STRUCT.iter_unpack(buf):
                    yield off, cds
                remaining -= len(buf) // _PAIR_STRUCT.size
        finally:
            f.close()

    return _gen()


def read_total_cds(source_path, *, filetype, compression):
    """Return the total CDS/protein count from the valid, fresh ``.didx`` header
    (O(1): a single header read, no body scan), or ``None`` to signal the caller to
    fall back to a full count. Used by ``domain_search -Z 0`` to learn the target
    count cheaply across (possibly sharded) databases."""
    validated = _validate_header(index_path_for(source_path), source_path,
                                 filetype=filetype, compression=compression)
    if validated is None:
        return None
    f, _n, total_cds = validated
    f.close()
    return total_cds


# --- shard naming / resolution -------------------------------------------------

def _split_db_name(db_path):
    """Split ``db_path`` into ``(dir, stem, fmt_ext, comp_ext)``.

    ``fmt_ext`` is the sequence-format extension (e.g. ``gb``/``fasta``); ``comp_ext``
    is a trailing compression extension or ``""``; ``stem`` is everything before
    ``fmt_ext`` (may itself contain dots). ``fmt_ext`` is ``""`` if there is no
    extension to anchor sharding on.
    """
    directory = os.path.dirname(db_path)
    base = os.path.basename(db_path)
    parts = base.split(".")
    comp = ""
    if len(parts) >= 2 and parts[-1].lower() in _COMPRESSION_EXTENSIONS:
        comp = parts[-1]
        parts = parts[:-1]
    if len(parts) >= 2:
        fmt = parts[-1]
        stem = ".".join(parts[:-1])
    else:
        fmt = ""
        stem = parts[0] if parts else ""
    return directory, stem, fmt, comp


def shard_path(db_path, shard_index=None, compress=None):
    """Build an output path derived from ``db_path``. ``shard_index`` inserts ``.N``
    after the base name (``None`` = no shard number, a single output). ``compress``
    overrides the output compression: ``True`` -> ``.bgz``, ``False`` -> none,
    ``None`` -> keep the input's compression suffix."""
    directory, stem, fmt, comp = _split_db_name(db_path)
    if compress is True:
        comp = "bgz"
    elif compress is False:
        comp = ""
    if shard_index is None:
        name = f"{stem}.{fmt}" if fmt else stem
    else:
        name = f"{stem}.{int(shard_index)}.{fmt}" if fmt else f"{stem}.{int(shard_index)}"
    if comp:
        name = f"{name}.{comp}"
    return os.path.join(directory, name) if directory else name


def shard_regex_for(db_path):
    """A regex matching the shard files of ``db_path`` in its directory. Capture
    group 1 is the shard number. Matches any compression suffix per shard (mixed
    compression across shards is allowed); does not match the unsharded file."""
    _directory, stem, fmt, _comp = _split_db_name(db_path)
    comp_alt = "|".join(re.escape(e) for e in sorted(_COMPRESSION_EXTENSIONS))
    fmt_part = re.escape(fmt)
    return re.compile(rf"^{re.escape(stem)}\.(\d+)\.{fmt_part}(?:\.(?:{comp_alt}))?$")


def resolve_database_shards(db_path):
    """Expand a logical database path to the ordered list of files to actually read.

    - only the unsharded file exists -> ``[db_path]``
    - only shards exist -> shard paths sorted numerically (leading zeros ignored)
    - both exist -> ``ValueError`` (ambiguous; refuse to guess)
    - neither exists -> ``[db_path]`` (let the normal open fail downstream)
    A leading-zero collision (``mydb.0.gb`` and ``mydb.00.gb``) raises; a gap in the
    shard numbering warns but proceeds.
    """
    db_path = str(db_path)
    base_exists = os.path.exists(db_path)
    directory, _stem, fmt, _comp = _split_db_name(db_path)

    shards_by_num = {}
    if fmt:
        pattern = shard_regex_for(db_path)
        search_dir = directory if directory else "."
        try:
            entries = os.listdir(search_dir)
        except OSError:
            entries = []
        for name in entries:
            m = pattern.match(name)
            if not m:
                continue
            num = int(m.group(1))
            full = os.path.join(directory, name) if directory else name
            if num in shards_by_num:
                raise ValueError(
                    f"Ambiguous shards for '{db_path}': {shards_by_num[num]} and {full} "
                    f"both resolve to shard number {num} (leading-zero collision)."
                )
            shards_by_num[num] = full

    if base_exists and shards_by_num:
        example = shards_by_num[min(shards_by_num)]
        raise ValueError(
            f"Both an unsharded database '{db_path}' and shard(s) (e.g. '{example}') exist; "
            "refusing to guess which to use. Remove one set."
        )

    if shards_by_num:
        nums = sorted(shards_by_num)
        if nums[0] != 0 or nums[-1] != len(nums) - 1:
            warnings.warn(
                f"Database shards for '{db_path}' are not a contiguous 0..N sequence "
                f"(found {nums}); proceeding with the shards present.",
                RuntimeWarning, stacklevel=2,
            )
        return [shards_by_num[n] for n in nums]

    return [db_path]
