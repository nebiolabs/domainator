"""Format domain_search databases: shard, BGZF-compress, and/or add `.didx` indexes.

A *domain_search database* is a GenBank or FASTA file, optionally BGZF-compressed,
optionally split into numbered shards. This tool prepares such databases:

- ``--shards N`` splits each input into N files (``mydb.gb`` -> ``mydb.0.gb`` ...),
  balanced by record count.
- ``--compress`` writes BGZF output (block-gzip: gunzip-readable, ~3-4x smaller, and
  still randomly seekable / partitionable), using the native parallel compressor.
- ``--index`` writes a ``.didx`` offset index for each output file so ``domain_search``
  skips the per-run offset scan.

At least one of ``--shards``, ``--compress``, or ``--index`` must be requested.
``--name`` renames the output database (keeping its format/compression extensions),
which is most useful alongside ``--shards``/``--compress``.
Operations compose (e.g. a plain ``.gb`` -> N BGZF shards each with a ``.didx``).
Shard content is produced by byte-exact range copy of the (decompressed) record
stream, so records are preserved exactly. Compressed inputs are decompressed to a
temporary working file first; uncompressed inputs are read in place.
"""
import os
import shutil
import sys
import warnings
from collections import namedtuple

from jsonargparse import ArgumentParser, ActionConfigFile

from domainator import __version__, RawAndDefaultsFormatter
from domainator import utils, db_index, bgzf_compress
from domainator.output_guardrails import (
    add_max_output_gb_argument,
    enforce_output_limit,
    make_temporary_output_path,
    max_output_gb_to_bytes,
)

_COPY_CHUNK = 1 << 20

# A single unit of work: produce `final_path` (and its index) from `[byte_start,
# byte_end)` of the uncompressed `working_path`. Picklable for the spawn pool.
ShardTask = namedtuple(
    "ShardTask",
    "working_path byte_start byte_end final_path out_bgzf needs_write "
    "do_index filetype level threads force",
)


def _copy_byte_range(src, dst, start, end):
    """Copy bytes ``[start, end)`` of ``src`` to a new file ``dst`` (streamed)."""
    with open(src, "rb") as fh, open(dst, "wb") as out:
        fh.seek(start)
        remaining = end - start
        while remaining > 0:
            chunk = fh.read(min(remaining, _COPY_CHUNK))
            if not chunk:
                break
            out.write(chunk)
            remaining -= len(chunk)


def _decompress_to_temp(input_path, tmp_path):
    """Decompress a BGZF/gzip input to an uncompressed temp file (streamed)."""
    handle, input_type = utils.open_seqfile(input_path, binary=True)
    try:
        with open(tmp_path, "wb") as out:
            while True:
                chunk = handle.read(_COPY_CHUNK)
                if not chunk:
                    break
                out.write(chunk)
    finally:
        if input_type == "name":
            handle.close()


def _record_boundaries(num_records, n_shards):
    """Split ``num_records`` into at most ``n_shards`` contiguous groups balanced by
    record count. Returns a list of ``(start_index, count)`` (length == effective
    number of shards)."""
    n_shards = max(1, min(n_shards, num_records))
    base, extra = divmod(num_records, n_shards)
    out = []
    start = 0
    for s in range(n_shards):
        count = base + (1 if s < extra else 0)
        out.append((start, count))
        start += count
    return out


def _format_shard_task(task):
    """Worker: write one final file (+ index) for a shard or whole-file output.
    Top-level so it is picklable under the spawn start method."""
    index_path = db_index.index_path_for(task.final_path)
    if not task.force and os.path.exists(task.final_path):
        if not task.do_index or os.path.exists(index_path):
            return task.final_path, "skipped"

    if task.needs_write:
        tmp = make_temporary_output_path(task.final_path)
        plain_tmp = None
        try:
            if task.out_bgzf:
                plain_tmp = tmp + ".plain"
                _copy_byte_range(task.working_path, plain_tmp, task.byte_start, task.byte_end)
                bgzf_compress.compress_to_bgzf(plain_tmp, tmp, level=task.level, threads=task.threads)
            else:
                _copy_byte_range(task.working_path, tmp, task.byte_start, task.byte_end)
            os.replace(tmp, task.final_path)
        finally:
            if os.path.exists(tmp):
                os.remove(tmp)
            if plain_tmp and os.path.exists(plain_tmp):
                os.remove(plain_tmp)

    if task.do_index:
        offsets, num_proteins = utils.get_offsets(task.final_path)
        compression = "bgzf" if task.out_bgzf else None
        db_index.write_index(task.final_path, offsets, num_proteins,
                             filetype=task.filetype, compression=compression)
    return task.final_path, ("written" if task.needs_write else "indexed")


def _plan_input(input_path, *, shards, compress, index, output_dir, name, level, force, log):
    """Plan the output tasks for one input. Returns ``(tasks, working_temp_or_None)``.
    Performs the per-input serial work (decompress + offset scan) up front so the
    returned tasks can be parallelized."""
    filetype = utils.get_file_type(input_path)
    if filetype not in ("genbank", "fasta"):
        raise ValueError(f"Input '{input_path}' is not a GenBank or FASTA file (got type {filetype!r}).")
    in_comp = utils.detect_compression(input_path)
    out_dir = output_dir or (os.path.dirname(input_path) or ".")
    os.makedirs(out_dir, exist_ok=True)
    if name is not None:
        # Rename the output database: keep the input's format + compression
        # extensions, swap the stem to `name`.
        _d, _stem, fmt, comp = db_index._split_db_name(input_path)
        new_basename = f"{name}.{fmt}" if fmt else name
        if comp:
            new_basename = f"{new_basename}.{comp}"
        base = os.path.join(out_dir, new_basename)
    else:
        base = os.path.join(out_dir, os.path.basename(input_path))

    # Preserve BGZF; uncompressed/gzip become BGZF only with --compress.
    out_bgzf = bool(compress) or (in_comp == "bgzf")

    # Fast path: a single output that is the input itself (no shard, no format or
    # location change) just gets an index in place -- no decompress, no scan, no copy.
    if shards is None:
        final = db_index.shard_path(base, shard_index=None, compress=out_bgzf)
        needs_write = (in_comp == "gzip") or (os.path.abspath(final) != os.path.abspath(input_path))
        if not needs_write:
            return [ShardTask(input_path, 0, 0, final, out_bgzf, False, index, filetype, level, 0, force)], None

    # Otherwise we copy record bytes, so we need an uncompressed working file.
    # Uncompressed inputs are read in place; compressed inputs decompress to a temp.
    working_temp = None
    if in_comp is None:
        working = input_path
    else:
        working_temp = make_temporary_output_path(base) + ".plain"
        print(f"Decompressing {input_path} -> temporary working file", file=log)
        _decompress_to_temp(input_path, working_temp)
        working = working_temp
    working_size = os.path.getsize(working)

    tasks = []
    if shards is None:
        # Single output that differs from the input (format/location change).
        tasks.append(ShardTask(working, 0, working_size, final, out_bgzf, True,
                               index, filetype, level, 0, force))
    else:
        if shards < 1:
            raise ValueError("--shards must be >= 1")
        offsets, _num_proteins = utils.get_offsets(working)   # only needed for shard boundaries
        num_records = len(offsets)
        if num_records == 0:
            if working_temp and os.path.exists(working_temp):
                os.remove(working_temp)
            raise ValueError(f"No sequence records found in '{input_path}'.")

        def byte_range(start_idx, count):
            b_start = int(offsets[start_idx])
            end_idx = start_idx + count
            b_end = int(offsets[end_idx]) if end_idx < num_records else working_size
            return b_start, b_end

        boundaries = _record_boundaries(num_records, shards)
        if len(boundaries) < shards:
            warnings.warn(
                f"'{input_path}' has only {num_records} records; producing {len(boundaries)} "
                f"shard(s) instead of {shards}.", RuntimeWarning, stacklevel=2)
        for i, (start_idx, count) in enumerate(boundaries):
            b_start, b_end = byte_range(start_idx, count)
            final = db_index.shard_path(base, shard_index=i, compress=out_bgzf)
            tasks.append(ShardTask(working, b_start, b_end, final, out_bgzf, True,
                                   index, filetype, level, 0, force))
        # Only a concern when the shards share the input's name (so domain_search,
        # given the input path, would resolve to both it and its shards). A --name
        # that renames the database avoids the ambiguity.
        _d, in_stem, _f, _c = db_index._split_db_name(input_path)
        same_name = name is None or name == in_stem
        if same_name and out_dir == (os.path.dirname(input_path) or ".") and os.path.exists(input_path):
            warnings.warn(
                f"Wrote shards alongside the original '{input_path}'. Remove the original "
                "(or use --output_dir) so domain_search does not see both the unsharded "
                "file and its shards.", RuntimeWarning, stacklevel=2)
    return tasks, working_temp


def format_db(inputs, *, shards=None, compress=False, index=False, output_dir=None,
              name=None, cpu=1, level=6, max_output_bytes=None, force=False, log=sys.stderr):
    """Shard / compress / index domain_search databases. Returns the list of final
    output paths produced (or confirmed). At least one of ``shards``, ``compress``,
    or ``index`` must be requested. ``name`` renames the output database (its format
    and compression extensions are preserved); it requires a single input."""
    if shards is None and not compress and not index:
        raise ValueError("Nothing to do: request at least one of --shards, --compress, or --index.")
    if name is not None:
        if len(inputs) > 1:
            raise ValueError("--name cannot be used with multiple inputs (each would collide on the same name).")
        if shards is None and not compress:
            warnings.warn(
                "--name on an index-only run copies the input to the new name before indexing "
                "(an index alone cannot rename a database); --name is intended for --shards/--compress.",
                RuntimeWarning, stacklevel=2)
    all_tasks = []
    working_temps = []
    try:
        for input_path in inputs:
            tasks, working_temp = _plan_input(
                input_path, shards=shards, compress=compress, index=index,
                output_dir=output_dir, name=name, level=level, force=force, log=log)
            all_tasks.extend(tasks)
            if working_temp:
                working_temps.append(working_temp)

        # Conservative upper bound: total bytes to write (compression only shrinks).
        projected = sum(t.byte_end - t.byte_start for t in all_tasks if t.needs_write)
        enforce_output_limit(
            projected_bytes=projected, max_output_bytes=max_output_bytes,
            output_description=f"domainator_format_db output ({len(all_tasks)} file(s))",
            mitigation_options=["--compress", "--shards"],
        )

        # One Rust compressor thread per pooled worker avoids cpu^2 oversubscription;
        # a single output uses all cores for that one (possibly large) file.
        single = len(all_tasks) == 1
        tasks = [t._replace(threads=(cpu if single else 1)) for t in all_tasks]

        if cpu > 1 and len(tasks) > 1:
            with utils.make_pool(processes=min(cpu, len(tasks))) as pool:
                results = pool.map(_format_shard_task, tasks)
        else:
            results = [_format_shard_task(t) for t in tasks]
    finally:
        for tmp in working_temps:
            if os.path.exists(tmp):
                os.remove(tmp)

    for path, status in results:
        print(f"{status}\t{path}", file=log)
    return [path for path, _ in results]


def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)
    parser.add_argument('-i', '--input', nargs='+', required=True, type=str,
                        help="GenBank or FASTA database file(s) to format. May be plain, gzip, or BGZF.")
    parser.add_argument('--shards', type=int, default=None,
                        help="Split each input into this many shards (balanced by record count). "
                             "Shards are named '<base>.0.<ext>', '<base>.1.<ext>', ... Default: no sharding.")
    parser.add_argument('--compress', '--bgzf', action='store_true', dest='compress', default=False,
                        help="Write BGZF (block-gzip) output (gunzip-readable, seekable, ~3-4x smaller).")
    parser.add_argument('--index', action='store_true', default=False,
                        help="Write a .didx offset index for each output file (lets domain_search "
                             "skip the per-run offset scan).")
    parser.add_argument('--output_dir', default=None, type=str,
                        help="Directory to write outputs into. Default: alongside each input.")
    parser.add_argument('--name', default=None, type=str,
                        help="Base name for the output database (format/compression extensions are kept), "
                             "e.g. --name mydb turns 'foo.gb' into 'mydb.gb' / 'mydb.0.gb.bgz'. Requires a "
                             "single input; intended for --shards/--compress (warns on index-only runs).")
    parser.add_argument('--cpu', type=int, default=1,
                        help="Number of parallel workers (also threads for the BGZF compressor).")
    parser.add_argument('--level', type=int, default=6,
                        help="BGZF/deflate compression level (1-12).")
    parser.add_argument('--force', action='store_true', default=False,
                        help="Overwrite existing outputs / rebuild indexes instead of skipping.")
    add_max_output_gb_argument(parser)
    parser.add_argument('--log', default=None, type=str,
                        help="Log file. Default: stderr.")
    parser.add_argument('--config', action=ActionConfigFile)

    params = parser.parse_args(argv)

    if params.cpu < 1:
        raise ValueError("--cpu must be >= 1")
    log = sys.stderr if params.log is None else open(params.log, "w")
    try:
        format_db(
            params.input, shards=params.shards, compress=params.compress, index=params.index,
            output_dir=params.output_dir, name=params.name, cpu=params.cpu, level=params.level,
            max_output_bytes=max_output_gb_to_bytes(params.max_output_gb),
            force=params.force, log=log,
        )
    finally:
        if params.log is not None:
            log.close()


def _entrypoint():
    main(sys.argv[1:])


if __name__ == '__main__':
    main(sys.argv[1:])
