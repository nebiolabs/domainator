[CLI reference](README.md)
# domainator infernal_build

Build a covariance model (CM) from a multiple sequence alignment (experimental).

```
Build a covariance model (CM) from a multiple sequence alignment (MSA).

NOTE: EXPERIMENTAL!

Allows the user to specify the ACC, NAME, and DESC fields of the CM profile.
Some of these options are not available in the standard cmbuild tool from Infernal.
The input alignment is passed directly to cmbuild, so any aligned MSA format it
accepts can be used here. Formats without structure annotation, such as aligned
FASTA, generally require --noss.
```

Give at least one of these outputs: `--cm_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE | required | Path of input aligned MSA. Any format accepted by cmbuild can be used. Formats without structure annotation, such as aligned FASTA, generally require --noss. Use - to read standard input. |
| `-o`, `--cm_out` | FILE |  | Write the CM to this file. Use - for standard output. |
| `--name` | STRING | required | Name of the CM profile. |
| `--acc` | STRING |  | Accession of the CM profile. |
| `--desc` | STRING |  | Description of the CM profile. |
| `--input_format` | STRING |  | Optional alignment format to pass through to cmbuild, such as stockholm, selex, or afa. |
| `--cmbuild_path` | PROGRAM |  | Path to the cmbuild executable. |
| `--noss` |  |  | Ignore secondary-structure annotation and build a zero-basepair model. |
| `--hand` |  |  | Use Stockholm RF annotation to define consensus columns. |
| `--symfrac` | NUMBER |  | Consensus-column threshold to pass through to cmbuild. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
