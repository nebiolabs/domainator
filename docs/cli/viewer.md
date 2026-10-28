[CLI reference](README.md)
# domainator viewer

Write the Domainator Viewer, empty or opening a network, tree, or projection bundle.

```
Write the Domainator Viewer, empty or opening a viewer bundle.

The viewer is one self-contained HTML page that opens the network, tree and projection
bundles (.dvb) that build_ssn, build_tree and build_projection write, and the sessions it
saves itself. Without -i/--input the page opens bundles you choose in it; with one, it opens
that bundle on load.
```

Give at least one of these outputs: `--html_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE |  | A Domainator Viewer bundle (.dvb) for the page to open on load: from build_ssn, build_tree or build_projection --dvb_out, or the viewer's Save session. Older network bundles (.dsnv) work too. Without it, the page opens the bundles you choose in it. |
| `-o`, `--html_out` | FILE |  | Write the viewer page to this HTML file. |
| `--name` | STRING |  | The page's title. Defaults to the bundle's name. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
