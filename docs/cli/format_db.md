[CLI reference](README.md)
# domainator format_db

Shard, compress, and index a database for domain_search.

```
Format domain_search databases: shard, BGZF-compress, and/or add `.didx` indexes.

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
```

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | GenBank or FASTA database file(s) to format. May be plain, gzip, or BGZF. |
| `--shards` | INTEGER |  | Split each input into this many shards (balanced by record count). Shards are named '<base>.0.<ext>', '<base>.1.<ext>', ... Default: no sharding. |
| `--compress` |  |  | Write BGZF (block-gzip) output (gunzip-readable, seekable, ~3-4x smaller). |
| `--index` |  |  | Write a .didx offset index for each output file (lets domain_search skip the per-run offset scan). |
| `--output_dir` | DIR |  | Directory to write outputs into. Default: alongside each input. |
| `--name` | STRING |  | Base name for the output database (format/compression extensions are kept), e.g. --name mydb turns 'foo.gb' into 'mydb.gb' / 'mydb.0.gb.bgz'. Requires a single input; intended for --shards/--compress (warns on index-only runs). |
| `--cpu` | INTEGER | 0 | Number of parallel workers (also threads for the BGZF compressor). 0: one per physical core. |
| `--level` | INTEGER | 6 | BGZF/deflate compression level (1-12). |
| `--force` |  |  | Overwrite existing outputs / rebuild indexes instead of skipping. |
| `--log_out` | FILE |  | Write the log to this file. If not supplied, writes it to stderr. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
