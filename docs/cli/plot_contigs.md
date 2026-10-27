[CLI reference](README.md)
# domainator plot_contigs

Draw contigs and their annotations as an interactive HTML page.

```
Create graphical representations of the contigs and their annotations

It works best for relatively small datasets, of a few megabases or a few hundred contigs.
```

Give at least one of these outputs: `--html_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | names of input genbank files. Use - to read standard input. |
| `-o`, `--html_out` | FILE | required | Write the plot to this html file. Use - for standard output. |
| `--height` | INTEGER | 2400 | height of the contig plot, in pixels |
| `--width` | INTEGER | 1200 | width of the contig plot, in pixels |
| `--types` | STRING ... |  | Feature types to include in the plot. Default is 'CDS exon repeat_region ncRNA rRNA tRNA Domainator Domain_Search'. |
| `--feature_height` | INTEGER | 50 | height of the feature arrows, in pixels. |
| `--feature_font_size` | INTEGER | 26 | font size of the feature labels, in pixels. |
| `--metadata` | FILE ... |  | tab separated files of additional contig metadata. |
| `--color_by` | STRING |  | Color the contig lines in the plot based on this column of the metadata table. |
| `--color_table` | FILE |  | tab separated file with two columns and no header: a --color_by value and a hex color, for example: CCDB   cc0000. Without it, colors are assigned automatically. |
| `--palette_offset` | INTEGER | 0 | Skip this many colors at the start of the palette before assigning colors to --color_by values. For example, if a plot colored by one column used the first 5 colors, --palette_offset 5 colors another column without reusing them. Ignored with --color_table. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
