[CLI reference](README.md)
# domainator system_db_download

Download system models (PADLOC), or convert a local copy, for --system_models.

```
Download system models, or convert a local copy, into a directory that other commands read

Converts a foreign system model database into Domainator's own system model format, and writes
it to --output_dir:

    systems.yaml     the models (see docs/file_formats.md)
    padloc.hmm       every profile the models use, named after the collection's database
    LICENSE          the source's license
    conversion.tsv   every approximation the conversion made: model, feature, what was done

The directory can then be given to --system_models (find_systems, domainate, select_by_system
...). Its profile file can be searched with domainate -r, and the hits are labeled the same way
as with domainate --system_models:

    domainator system_db_download --source padloc --output_dir padloc_db
    domainator domainate -i genomes.gb -r padloc_db/padloc.hmm -e 1e-3 -o annotated.gb
    domainator find_systems -i annotated.gb --system_models padloc_db -o systems.gb

Sources:
    padloc: PADLOC's padloc-db (https://github.com/padlocbio/padloc-db, MIT license). Without
        --input, a release is downloaded from GitHub (--release, default: the latest). With
        --input, a local copy (a directory with sys/, hmm/, hmm_meta.txt and sys_meta.txt) is
        converted. CRISPR arrays can't be detected, so they are dropped from the models.
```

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `--source` | {padloc} | required | The format of the database to download or convert. |
| `-i`, `--input` | FILE |  | A local copy of the database to convert, instead of downloading it. |
| `--release` | STRING |  | The release to download, e.g. v2.0.0. Default: the latest, which is recorded in systems.yaml. Not with --input. |
| `--output_dir` | DIR | required | The directory to write the models, profiles, license and conversion report to. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
