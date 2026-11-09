[CLI reference](README.md)
# domainator import_systems

Convert system calls made by other tools (antiSMASH) into Domain_System features.

```
Convert system calls made by other tools into Domain_System features

Reads GenBank records that another tool has annotated, and adds a Domain_System feature (see
find_systems) for each system it called, so that Domainator's selection, reporting, coloring and
plotting tools work on them.

Sources:
    antismash: antiSMASH GenBank output. --level protocluster (the default) makes one system per
        protocluster, located at its core, named by its product, in the antismash database; its
        members are the core CDSs whose profiles triggered the product. --level candidate_cluster
        and --level region make one system per candidate cluster (with antiSMASH's /kind) or
        region, located at the whole feature. Systems imported before are replaced.

    domainator import_systems -i antismash_output.gbk --source antismash -o systems.gb
    domainator select_by_system -i systems.gb --system_categories PKS --model_neighborhood -o pks_regions.gb
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | GenBank files written by the source tool. Use - to read standard input. |
| `-o`, `--genbank_out` | FILE |  | Write the records to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the records to this FASTA file. Use - for standard output. |
| `--source` | {antismash} | antismash | The tool that made the system calls. |
| `--level` | {protocluster,candidate_cluster,region} | protocluster | antiSMASH only: which antiSMASH feature becomes a system. protocluster: one per protocluster, located at its core. candidate_cluster, region: one per candidate cluster or region, located at the whole feature. |
| `--clear_source` |  |  | Remove the source tool's cluster features (antiSMASH's region, cand_cluster, protocluster and proto_core) after converting them. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
