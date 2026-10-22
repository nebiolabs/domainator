[CLI reference](README.md)
# domainator color_table_to_legend

Draw an SVG legend from a color table.

```
Generates an SVG legend based on a color table

The color table is a tsv file that maps annotations to hex color codes.
The generated SVG file contains rectangles filled with the corresponding colors and text labels for each annotation.
```

Give at least one of these outputs: `--svg_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | names of color table files.                         Files are tab separated with two columns and no header, columns are: annotation, hex color. For example: CCDB   cc0000 Use - to read standard input. |
| `-o`, `--svg_out` | FILE | required | Write the legend to this svg file. Use - for standard output. |
| `--title` | STRING | Legend | Title of the legend. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
