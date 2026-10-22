"""Reading and writing GenBank and FASTA files.

Always read and write sequence files through :func:`parse_seqfiles` and :func:`write_genbank`,
which keep the qualifiers Domainator adds.
"""
import copy
from datetime import date
from pathlib import Path
import gzip
import os
import sys
import warnings
from typing import IO, Container, Iterable, Iterator, Optional, Union

from domainator import lean_record
from Bio import SeqIO
from Bio import bgzf
from Bio.SeqRecord import SeqRecord
from domainator.genbank import parse_genbank, write_genbank_records
from domainator.model import get_domain_cds_annotation_strings
from domainator.partial import normalize_partial_codes
try:
    from domainator import _gbfast  # optional native (Rust) acceleration
except ImportError:
    _gbfast = None


EXTENSION_TO_TYPE = {
    "gb":"genbank",
    "gbk":"genbank",
    "genbank":"genbank",
    "gbff":"genbank",
    "fasta":"fasta",
    "fa":"fasta",
    "fna":"fasta",
    "pep":"fasta",
    "faa":"fasta",
    # "gff":"gff",
    # "gtf":"gff",
    # "gff3":"gff",
    "hmmer":"hmm",
    "hmm":"hmm",
    "cm":"cm",
    "hdf5":"hdf5",
    "hdf":"hdf5",
    "h5":"hdf5",
    # Protein structure files, read by the structure_* tools via an external structural
    # aligner. No sequence reader here handles them: every consumer of get_file_type
    # whitelists the types it accepts, so these are recognized only so that callers can say
    # "that is a structure file" instead of "extension not recognized".
    "pdb":"structure",
    "pdb1":"structure",
    "ent":"structure",
    "cif":"structure",
    "mmcif":"structure",
    "bcif":"structure",
}

# Compression suffixes that are stripped to reveal the underlying sequence format
# (e.g. "foo.gb.gz" -> genbank). ".bgz"/".bgzf" indicate BGZF (random-access gzip).
COMPRESSION_EXTENSIONS = {"gz", "bgz", "bgzf"}

GZIP_MAGIC = b"\x1f\x8b"


STDIN_PATHS = frozenset({"-", "/dev/stdin"})
"""Input paths that mean "read standard input"."""
STDOUT_PATHS = frozenset({"-", "/dev/stdout"})
"""Output paths that mean "write standard output"."""


class _UnclosedStdout:
    """Standard output, wrapped so that a tool can close it like any output file: closing only flushes."""

    def __init__(self, stream):
        self._stream = stream

    def write(self, data):
        return self._stream.write(data)

    def flush(self):
        self._stream.flush()

    def close(self):
        self._stream.flush()

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()

    def __getattr__(self, name):
        return getattr(self._stream, name)

def db_name_from_path(name) -> str:
    """Derive the database name Domainator labels annotations with from a file path.

    Strips a trailing compression suffix before the format suffix, so that
    "refs.hmm" and "refs.hmm.gz" both yield "refs". Without this, a compressed
    reference would be labeled "refs.hmm" and produce output that differs from
    the uncompressed one.

    Only the format suffix is removed after the compression suffix, and only when
    something is left over, so paths with unknown-but-meaningful extensions keep
    behaving as they always have ("some.long.name.hmm" -> "some.long.name",
    "refs.v1" -> "refs", "db.gz" -> "db").
    """
    stem = os.path.basename(os.fspath(name))
    suffix = Path(stem).suffix
    if suffix and suffix[1:].lower() in COMPRESSION_EXTENSIONS:
        stem = stem[: -len(suffix)]
    trimmed = Path(stem).stem
    return trimmed if trimmed else stem


def get_file_type(filename):
    out = None
    suffixes = [s[1:].lower() for s in Path(filename).suffixes]
    # Strip a trailing compression suffix so "foo.gb.gz" is treated as genbank.
    if suffixes and suffixes[-1] in COMPRESSION_EXTENSIONS:
        suffixes = suffixes[:-1]
    if suffixes:
        file_extension = suffixes[-1]
        if file_extension in EXTENSION_TO_TYPE:
            out = EXTENSION_TO_TYPE[file_extension]
    return out


def detect_compression(path):
    """Inspect magic bytes to classify a file as 'bgzf', 'gzip', or None.

    BGZF blocks always start with b"\\x1f\\x8b\\x08\\x04" (gzip magic with the
    FEXTRA flag set); plain gzip has the gzip magic but a different flag byte.
    Returns None for uncompressed files (or if the path can't be read).
    """
    try:
        with open(path, "rb") as fh:
            magic = fh.read(4)
    except (OSError, TypeError, ValueError):
        return None
    if magic[:2] != GZIP_MAGIC:
        return None
    if magic == b"\x1f\x8b\x08\x04":
        return "bgzf"
    return "gzip"


def open_seqfile(filename_or_handle, binary=False):
    """Open a sequence file for reading, transparently decompressing gzip/BGZF.

    Returns (handle, input_type) where input_type is "name" (we opened it) or
    "handle" (an already-open object was passed through unchanged).

    BGZF inputs are opened via BgzfReader, whose tell()/seek() use 64-bit virtual
    offsets, so the offset-partition scheme keeps working for parallel random
    access. Plain (non-blocked) gzip is supported as a convenience but its seeks
    re-decompress from the start of the file, which defeats parallel seeking, so
    BGZF is the recommended compressed database format.
    """
    if not isinstance(filename_or_handle, (str, os.PathLike)):
        return filename_or_handle, "handle"

    compression = detect_compression(filename_or_handle)
    if compression == "bgzf":
        mode = "rb" if binary else "r"
        return bgzf.BgzfReader(str(filename_or_handle), mode=mode), "name"
    elif compression == "gzip":
        mode = "rb" if binary else "rt"
        return gzip.open(filename_or_handle, mode), "name"
    else:
        mode = "rb" if binary else "r"
        return open(filename_or_handle, mode), "name"

def is_compressed_path(path):
    """Return True if the path's extension marks it as gzip/BGZF (.gz/.bgz/.bgzf)."""
    if not isinstance(path, (str, os.PathLike)):
        return False
    suffixes = [s[1:].lower() for s in Path(path).suffixes]
    return bool(suffixes) and suffixes[-1] in COMPRESSION_EXTENSIONS


def open_writable_seqfile(path, compressed=None):
    """Open a path for writing, using a BGZF writer for compressed outputs.

    ``-`` and ``/dev/stdout`` return standard output, which is never compressed and which
    ``close()`` only flushes.

    BGZF output is itself valid gzip (so ``gunzip`` reads it) but is also
    randomly seekable, so the file can be re-used directly as a compressed,
    partitionable input to the next domainator step. ``compressed`` overrides the
    extension-based decision (useful when writing to a temp path whose own
    extension does not reflect the final destination).
    """
    if path in STDOUT_PATHS:
        return _UnclosedStdout(sys.stdout)
    if compressed is None:
        compressed = is_compressed_path(path)
    if compressed:
        return bgzf.BgzfWriter(str(path))
    return open(path, "w")


def open_if_is_name_for_write(filename_or_handle):
    """Like open_if_is_name(mode="w") but routes compressed paths through BGZF.

    Returns (handle, input_type) where input_type is "name" (we opened it) or
    "handle" (an already-open object passed through unchanged).
    """
    if not isinstance(filename_or_handle, (str, os.PathLike)):
        return filename_or_handle, "handle"
    if filename_or_handle in STDOUT_PATHS:
        return sys.stdout, "handle"
    if is_compressed_path(filename_or_handle):
        return bgzf.BgzfWriter(str(filename_or_handle)), "name"
    return open(filename_or_handle, "w"), "name"


def open_if_is_name(filename_or_handle, mode="r"):
    """
        if a file handle is passed, return the file handle
        if a Path object or path string is passed, open and return a file handle to the file.

        returns:
            file_handle, input_type ("name" | "handle")
    """
    out = filename_or_handle
    input_type = "handle"
    try:
        out = open(filename_or_handle, mode)
        input_type = "name"
    except TypeError:
        pass
    except Exception as e:
        raise(e)

    return (out, input_type)


def write_genbank(seq_records, file_name, default_molecule_type="DNA"):
    """Write SeqRecords as GenBank, to a file name or a handle. The records are not changed.

    Each record is written as parse_seqfiles would read it back: its id goes in the LOCUS line
    (with spaces removed) and its name in ACCESSION. Records without a molecule type get
    default_molecule_type, and records without a date get today's. Features are written in order
    of start position, and each CDS gets a domainator_<database> qualifier summarizing the
    Domainator domains on it (and loses any for databases with no domains left).
    """
    outfile, file_type = open_if_is_name_for_write(file_name)
    try:
        for record in seq_records:
            write_genbank_records([_genbank_output_record(record, default_molecule_type)], outfile)
    finally:
        if file_type == "name":
            outfile.close()


def _genbank_output_record(record, default_molecule_type):
    """A new record that is ``record`` as write_genbank writes it, sharing what it doesn't change."""
    annotations = dict(record.annotations)
    annotations.setdefault("molecule_type", default_molecule_type)
    if "date" not in annotations:
        annotations["date"] = date.today().strftime("%d-%b-%Y").upper()
    cds_domain_annotations = get_domain_cds_annotation_strings(record)
    features = [_with_domain_qualifiers(f, cds_domain_annotations) if f.type == "CDS" else f for f in record.features]
    features.sort(key=lambda x: x.location.start)
    return SeqRecord(record.seq, id=record.name, name=record.id.replace(" ", ""), description=record.description,
                     dbxrefs=record.dbxrefs, features=features, annotations=annotations,
                     letter_annotations=record.letter_annotations)


def _with_domain_qualifiers(cds, cds_domain_annotations):
    """``cds``, or a copy of it with its domainator_<database> qualifiers brought up to date."""
    domains = cds_domain_annotations.get(cds.qualifiers.get("cds_id", [None])[0], {})
    wanted = {f"domainator_{db}": [" & ".join(strings)] for db, strings in domains.items()}
    stale = [q for q in cds.qualifiers if q.startswith("domainator_") and q not in wanted]
    if not stale and all(cds.qualifiers.get(q) == v for q, v in wanted.items()):
        return cds
    qualifiers = {k: v for k, v in cds.qualifiers.items() if k not in stale}
    qualifiers.update(wanted)
    out = copy.copy(cds)
    out.qualifiers = qualifiers
    return out


def swap_name_id(record):
    """
        necesarry because biopython handles genbank names and ids differently than other programs, such as Geneious
    """
    id = record.id
    record.id = record.name
    record.name = id

def copy_SeqRecord(seqrecord):
    new_record = SeqRecord(seq=seqrecord.seq, id=seqrecord.id, name=seqrecord.name,
                           description=seqrecord.description,
                           dbxrefs=seqrecord.dbxrefs.copy(),
                           features=seqrecord.features.copy(),
                           annotations=seqrecord.annotations.copy(),
                           letter_annotations=seqrecord.letter_annotations,
                           )
    return new_record


def native_parser_available():
    """True when the optional native (Rust) ``_gbfast`` parser/scanner is installed."""
    return _gbfast is not None


def _resolve_genbank_parser(requested):
    """Decide which GenBank parser backend to use ("lean" or "biopython").

    Two tiers: "lean" uses the optional native (Rust) ``_gbfast`` extension and yields
    lean record objects, so only callers that handle those request it; "biopython" uses
    the Biopython parser (no build dependency). "lean" downgrades to "biopython"
    when the native extension is absent.
    """
    backend = (requested or "biopython").lower()
    if backend not in ("biopython", "lean"):
        raise ValueError(f"genbank_parser must be 'biopython' or 'lean', not {requested!r}.")
    if backend == "lean" and _gbfast is None:
        backend = "biopython"
    return backend


# Paths for which we've already warned that the native parser fell back to
# Biopython. The fallback is expected for some inputs (protein GenBank,
# nonstandard LOCUS lines, '*'-containing sequences), and the per-partition parse
# would otherwise emit one warning per partition, so warn at most once per file
# (per process; spawned workers each warn once, which is still bounded).
_native_fallback_warned_paths = set()


def _warn_native_fallback(file):
    key = str(file)
    if key in _native_fallback_warned_paths:
        return
    _native_fallback_warned_paths.add(key)
    warnings.warn(
        f"The native (Rust) parser could not fully parse '{file}'; falling back to "
        "the Biopython parser for this file. This is expected for protein GenBank, "
        "nonstandard LOCUS lines, and '*'-containing sequences; results are "
        "unaffected (only this file's parsing is slower).",
        RuntimeWarning,
        stacklevel=2,
    )


def _native_usable(file):
    """True when the native (Rust) parser can be used for this file.

    Applies to uncompressed real paths and to BGZF (the Rust parser seeks by the
    64-bit virtual offset the Python offset scanner produces). Plain (non-blocked)
    gzip has no random access, so it stays on the Biopython path.
    """
    if _gbfast is None or not isinstance(file, (str, os.PathLike)):
        return False
    return detect_compression(file) in (None, "bgzf")


def _iter_raw_genbank_records(file, seek_to, default_molecule_type, backend, max_recs=float("inf")):
    """Yield pre-(name/id)-swap records for one GenBank file.

    "lean" uses the native (Rust) parser to yield LeanContig/LeanSearchContig
    objects and transparently falls back to the Biopython parser (resuming after
    the records the native parser already emitted) for records it cannot handle
    (e.g. a nonstandard LOCUS line or a sequence it rejects). It only applies to
    uncompressed real paths; compressed/handle inputs use Biopython. Handles are
    closed in a finally so abandoning the iterator early (e.g. on a max_recs
    break) does not leak file descriptors.
    """
    if backend == "lean" and _native_usable(file):
        compressed = detect_compression(file) == "bgzf"

        def produce():
            bound = int(max_recs) if max_recs != float("inf") else -1
            yield from lean_record.iter_lean_search_native(
                file, seek_to or 0, bound, default_molecule_type, compressed
            )

        emitted = 0
        try:
            for rec in produce():
                emitted += 1
                yield rec
            return
        except lean_record.LeanParseError:
            _warn_native_fallback(file)
        # Native parser stopped early: re-read with Biopython, converting to lean
        # and skipping the records already emitted.
        bp_handle, bp_input_type = open_seqfile(file)
        try:
            if seek_to is not None:
                bp_handle.seek(seek_to)
            for skip_count, rec in enumerate(parse_genbank(bp_handle)):
                if skip_count < emitted:
                    continue
                yield lean_record.seqrecord_to_lean(rec)
        finally:
            if bp_input_type == "name":
                bp_handle.close()
    else:
        handle, input_type = open_seqfile(file)
        try:
            if seek_to is not None:
                handle.seek(seek_to)
            for rec in parse_genbank(handle):
                yield rec
        finally:
            if input_type == "name":
                handle.close()


def _iter_raw_fasta_records(file, seek_to, backend="biopython", default_molecule_type=None, max_recs=float("inf")):
    """Yield records for one FASTA file.

    With backend == "lean" and the native parser available, yields LeanFastaContig
    objects (the bulk fast path). Otherwise yields Biopython SeqRecords.
    """
    if backend == "lean" and _native_usable(file) and getattr(lean_record, "LeanFastaContig", None) is not None:
        compressed = detect_compression(file) == "bgzf"
        bound = int(max_recs) if max_recs != float("inf") else -1
        yield from lean_record.iter_lean_fasta_native(file, seek_to or 0, bound, default_molecule_type, compressed)
        return
    handle, input_type = open_seqfile(file)
    try:
        if seek_to is not None:
            handle.seek(seek_to)
        for rec in SeqIO.parse(handle, "fasta"):
            yield rec
    finally:
        if input_type == "name":
            handle.close()


def parse_seqfiles(seqfiles: Iterable[Union[str, os.PathLike, IO]], contigs: Optional[Container[str]] = None,
                   filetype_override: Optional[str] = None, seek_to: Optional[int] = None,
                   max_recs: float = float("inf"), default_molecule_type: str = 'protein',
                   genbank_parser: Optional[str] = None) -> Iterator[SeqRecord]:
    """Read the records of GenBank and FASTA files, one at a time.

    Args:
        seqfiles: GenBank or FASTA paths or open handles. "-" or "/dev/stdin" reads standard
            input, as GenBank unless ``filetype_override`` says otherwise.
        contigs: yield only the records with these ids (any container, such as a set or the
            keys of a dict). None yields every record.
        filetype_override: "genbank" or "fasta", instead of guessing the format from the file
            extension. Required for handles.
        seek_to: seek to this byte offset in the file before reading any records.
        max_recs: stop reading each file after yielding this many of its records.
        default_molecule_type: the molecule_type of records that don't declare one, such as
            FASTA records.
        genbank_parser: which GenBank backend to use, "lean" or "biopython" (the default).
            "lean" uses the optional native (Rust) parser and yields LeanContig objects; it is
            much faster, but falls back to Biopython for records it cannot parse and when the
            native extension is not built.

    Yields:
        The records, as SeqRecords (or LeanContigs with ``genbank_parser="lean"``).
    """


    for file in seqfiles:
        file_type = filetype_override
        if file in STDIN_PATHS:
            file = sys.stdin
            file_type = filetype_override or "genbank"
        is_handle = not isinstance(file, (str, os.PathLike))

        if file_type is None:
            if is_handle:
                raise ValueError(f"filetype_override must be specified when passing a handle.")
            file_type = get_file_type(file)
            if file_type == "structure":
                raise ValueError(
                    f"'{file}' is a protein structure file, not a sequence file. Convert it "
                    "with domainator structure_to_genbank first, or use domainator structure_domainate / "
                    "domainator structure_search to work with structures directly."
                )
            if file_type not in {"genbank", "fasta"}:
                raise ValueError(f"file extension not recognized: {file}, please specify format.")

        # the native parser needs a real, seekable binary path; arbitrary handles
        # stay on the Biopython parser.
        backend = "biopython" if is_handle else _resolve_genbank_parser(genbank_parser)
        if file_type == "genbank":
            raw_iter = _iter_raw_genbank_records(file, seek_to, default_molecule_type, backend, max_recs=max_recs)
        else:
            raw_iter = _iter_raw_fasta_records(file, seek_to, backend, default_molecule_type, max_recs)

        lean_search_types = lean_record.LEAN_SEARCH_TYPES if lean_record is not None else ()
        recs_returned = 0
        try:
            for rec in raw_iter:
                if recs_returned >= max_recs:
                    break
                # Native search records (GenBank/FASTA) are already post-swap and
                # carry molecule_type internally, so they skip the swap /
                # molecule_type fallback below.
                is_search_contig = bool(lean_search_types) and isinstance(rec, lean_search_types)
                if file_type == "genbank" and not is_search_contig:
                    swap_name_id(rec)
                    normalize_partial_codes(rec)
                if contigs is None or rec.id in contigs:
                    if not is_search_contig and 'molecule_type' not in rec.annotations:  # TODO: is this ok (will a nucleotide sequence always have molecule_type defined and a peptide sequence always not have?)
                        rec.annotations['molecule_type'] = default_molecule_type
                    yield rec
                    recs_returned += 1
        finally:
            # Closing the sub-generator runs its finally block, closing the handle
            # even when we broke out early on max_recs.
            raw_iter.close()
