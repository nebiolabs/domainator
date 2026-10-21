[index](README.md)
# The tool manifest (`domainator describe`)

`domainator describe` prints one JSON object that describes every command: its options, what
kind of value each one takes, and which ones are files to read or write. It is built from the
commands' own parsers, so it always matches the command line. It is meant for programs that
drive Domainator, such as a GUI or a workflow system, which render a form for each command from
it instead of keeping their own copy of the options.

```
domainator describe --indent 1 > manifest.json
domainator describe --tool domainate select_by_cds
```

## Running a command from the manifest

Write the chosen values to a JSON file whose keys are the parameters' `name`s, then run the
command with it:

```
domainator domainate --config options.json
```

```json
{"input": ["genomes.gb"], "reference": ["Pfam-A.hmm"], "genbank_out": "annotated.gb", "cpu": 8}
```

Leave out a parameter to get its default. `-` in an input or output path means standard input
or output, where the parameter's `stdin` or `stdout` is `true`.

## Schema (manifest version 1)

```javascript
{
  "manifest_version": 1,           // see "Versions" below
  "domainator_version": "...",
  "invocation": {"command": ["domainator", "{tool}"], "config_option": "--config",
                 "config_keys": "parameter name", "stdio_path": "-"},
  "categories": [{"name": "annotate", "title": "..."}, ...],   // in display order
  "tools": [
    {
      "name": "domainate",              // the command: domainator domainate
      "category": "annotate",
      "summary": "...",                 // one line
      "description": "...",             // the command's documentation
      "parameters": [ ... ],            // in the order of --help
      "required_outputs": ["genbank_out", "fasta_out"],   // at least one of these must be given; null if none is required
      "mutually_exclusive": [["verbose", "quiet"], ...]   // at most one of each group
    }
  ]
}
```

Each parameter:

| Key | Meaning |
|---|---|
| `name` | The config key (the option's name without dashes, e.g. `max_hits_per_cds`). |
| `flags` | The command-line flags, e.g. `["-e", "--evalue"]`. |
| `help` | The option's help text. |
| `kind` | `flag` (on or off), `value`, `path`, or `columns` (see below). |
| `type` | `string`, `integer`, `number` or `boolean`. |
| `min`, `max`, `exclusive_min` | Bounds, when the option has them. |
| `role` | For a `path`: `input` (a file to read), `output` (a file to write), `directory`, or `program` (an external program's executable). |
| `multiple` | `true` if the option takes a list of values. |
| `count` | For an option that takes exactly this many values (more than one). |
| `required` | `true` if the command needs it. |
| `default` | The default value. Missing when `default_dynamic` is `true`: the default depends on the user's environment (e.g. a directory in their home), so leave the option out to get it. |
| `choices` | The allowed values, when there is a fixed list. |
| `stdin` | For `input`: whether `-` (standard input) is accepted. Commands that read their input more than once need a file. |
| `stdout` | For an output: whether `-` (standard output) is accepted. At most one output can be standard output. |
| `sets` | For a `flag`: the value it sets when given (`true` for `--x`, `false` for `--no_x`). |
| `negative_flags` | For a switch, the flags that turn it off (`--no_<name>`). |

A `columns` parameter (the `--columns` option of `enum_report` and `hmmer_report`) is an ordered
list of column tokens. `columns` lists the available columns: each has a `name`, the names of its
`args` (written after the name, separated by `:`, as in `qualifier:CDS:product`), `repeat_last`
(the last argument can be repeated, as in `taxid:genus:species`), and `help`.

How options are classified: an option is an output exactly when its name ends in `_out`. Names
ending in `_file`, and `input`, `reference`, `metadata` and a few others, are inputs; names ending
in `_dir` are directories; names ending in `_path` are external programs (except
`ncbi_taxonomy_path`, a directory). `test/test_manifest.py` checks that every option of every
command is classified.

## Versions

`manifest_version` is incremented when a key is removed or changes meaning. New keys can
appear without a new version, so readers must ignore keys they don't know. The commands' options
themselves follow the 1.x stability contract ([stability.md](stability.md)): within 1.x an
option keeps its name and meaning, and new options can appear.

## Notes for GUI authors

Domainator 0.x bundled a browser GUI (`domainator_server`), which is being rebuilt as a
separate project on top of this manifest. A review of the old server found problems that a new
GUI that runs Domainator for other users should avoid:

- It passed unknown parameters straight through to the command line, which, together with the
  `program` parameters (paths to executables), allowed running arbitrary programs. Accept only
  the parameters the manifest lists, and don't expose `role: program` parameters to untrusted
  users.
- It accepted any path on the server as a file parameter. Restrict `role: input` paths to files
  the user uploaded, and choose output paths on the server side.
- It served uploaded HTML from its own origin, so an uploaded page could run scripts in the GUI
  (stored cross-site scripting). Serve user files as downloads, or from a separate origin.
- It had no authentication and no Host/Origin check, its `--debug` mode turned on the Werkzeug
  debugger (remote code execution), and it didn't validate file ids.

Policy that belongs in the GUI rather than in Domainator: which commands and parameters to show
(and which to tuck away as advanced), presets that differ from the command-line defaults,
resource limits, and sandboxing.
