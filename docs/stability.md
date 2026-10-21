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

The on-disk formats in [file_formats.md](file_formats.md): the GenBank feature types and
qualifiers Domainator writes, matrix hdf5 files, `.dsnv` bundles and `.didx` indexes. Each
versioned format is read for its current and previous major versions, and a file from a newer
major version is refused with a `FormatVersionError`.
