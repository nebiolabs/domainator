[CLI reference](README.md)
# domainator system_db_download

Download system models (PADLOC), or convert a local copy, for --system_models.

```
Download system models, or convert a local copy, into a directory that other commands read

Converts a foreign system model database into Domainator's own system model format, and writes
it to --output_dir:

    systems.yaml     the models (see docs/file_formats.md)
    <database>.hmm   every profile the models use, named after the collection's database
                     (padloc, defensefinder, or a MacSyFinder package's name)
    LICENSE          the source's license
    conversion.tsv   every approximation the conversion made: model, feature, what was done

The directory can then be given to --system_models (find_systems, domainate, system_search
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
    defensefinder: DefenseFinder's defense-finder-models (https://github.com/mdmparis/defense-finder-models;
        its metadata gives the licenses GPL-3.0 and CC BY-NC-SA 4.0). Without --input, a release
        is downloaded from GitHub. With --input, a local copy (a MacSyFinder package directory,
        with metadata.yml, definitions/ and profiles/, as installed by defense-finder update or
        macsydata) is converted. Hits must cover --coverage_profile of a profile (default 0.4,
        as DefenseFinder uses) and pass its gathering threshold.
    macsyfinder: any other MacSyFinder model package (models format 2.0), given with --input. It
        is named after its directory. --coverage_profile defaults to the package's
        model_conf.xml, else MacSyFinder's 0.5.

MacSyFinder loners and multi-locus models are found as single loci, and listed in
conversion.tsv.
```

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `--source` | {padloc,defensefinder,macsyfinder} | required | The format of the database to download or convert. |
| `-i`, `--input` | FILE |  | A local copy of the database to convert, instead of downloading it. |
| `--release` | STRING |  | The release to download, e.g. v2.0.0 (padloc) or 3.1.0 (defensefinder). Default: the latest, which is recorded in systems.yaml. Not with --input. |
| `--coverage_profile` | NUMBER |  | For defensefinder and macsyfinder: the fraction of a profile a hit must cover. Default: 0.4 for defensefinder (as DefenseFinder runs MacSyFinder); for macsyfinder, the package's model_conf.xml, else 0.5. |
| `--output_dir` | DIR | required | The directory to write the models, profiles, license and conversion report to. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
