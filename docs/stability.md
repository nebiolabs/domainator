# What Domainator 1.x keeps stable

From 1.0 on, a change that breaks any of the following waits for 2.0. Anything not listed is
internal and can change in any release.

## Commands

Every `domainator <command>`, its options (names, meanings and defaults), and its exit
status: 0 on success, 1 for an error the user can fix (printed as one line), 2 for a usage
error.

## Python API

- The names the package root exports, listed in `domainator.__all__`: `parse_seqfiles`,
  `write_genbank`, `DataMatrix`, the exceptions in `domainator.errors`, the feature-type
  constants `DOMAIN_FEATURE_NAME` and `DOMAIN_SEARCH_BEST_HIT_NAME`, and `__version__`.
- Each command's function, `from domainator.<command> import <command>`. Its first parameter
  is positional and the rest are keyword-only, named like the command's options, so new
  parameters can be added without breaking callers. What it returns (records, a `DataMatrix`,
  ...) is part of the contract.
- The other modules (`seqio`, `hmmio`, `model`, `locations`, `partial`, `taxonomy`, `expr`,
  `colors`, `offsets`, `external`, `cli`, ...) and every name that starts with `_` are internal.

## File formats

The formats in [file_formats.md](file_formats.md):

- **GenBank**: the `Domainator` and `Domain_Search` feature types and the qualifiers in its
  qualifier reference, including how locations and coordinates are written. These change only
  by addition, so a reader must ignore qualifiers it doesn't know.
- **Versioned files**: matrix hdf5 files (`MATRIX_FILE_VERSION`), `.dsnv` bundles and `.didx`
  indexes. Each is read for its current and previous major versions, and a file from a newer
  major version is refused with a `FormatVersionError`.
- **Report JSON** (`--json_out` of `summary_report`, `matrix_report` and `ssn_navigator`): the
  keys, versioned by the `schema_version` key. New keys can appear in 1.x.
- **The tool manifest** (`domainator describe`, [tool_manifest.md](tool_manifest.md)), versioned
  by its `manifest_version` key. New keys can appear in 1.x.
- **Other outputs**: the columns of the report TSVs and NDJSON, `hmmer_compare` scores and
  alignments, XGMML attributes and Newick naming, and the color-table and dense-text matrix
  formats.

What is not part of the contract: HTML reports and the SSN viewer page (their look and the
data embedded in them), text summaries meant for reading, the order of hits that tie on score,
warning and log messages, and the floating-point details of scores (a different HMMER, DIAMOND
or Foldseek version can change them).
