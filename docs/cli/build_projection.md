[CLI reference](README.md)
# domainator build_projection

Project a matrix into 2D or 3D with UMAP, t-SNE, or PCA.

```
Generates 2d or 3d umap, tsne, or pca projections from matrices

    output can be images, coordinates, or Cytoscape xgmml files.
```

Give at least one of these outputs: `--xgmml_out`, `--png_out`, `--coords_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE | required | The matrix that you want to project into 2D or 3D space, rows are observations (typically gene names), columns are features (typically similarity to other genes, or hmm profile scores). Format can be tab-separated text or Domainator hdf5. |
| `--algorithm` | {pca,tsne,umap,coords} | pca | Which projection algorithm to use. |
| `--metadata` | FILE ... |  | tab separated files of sequence metadata. |
| `--color_by` | STRING |  | Color the points in the output image based on this column of the metadata table. |
| `--color_table` | FILE |  | tab separated file with two columns and no header, columns are: annotation, hex color. For example: CCDB   cc0000 |
| `--color_table_out` | FILE |  | tab separated file with two columns and no header, columns are: annotation, hex color. Written after the color table is updated with new colors, for example if using --color_by, but not supplying an external color table. |
| `-o`, `--xgmml_out` | FILE |  | Write the projection to this Cytoscape xgmml file. |
| `--png_out` | FILE |  | Write a picture of the projection to this png file. |
| `--coords_out` | FILE |  | Write the projected coordinates to this tab-separated file. |
| `--pca_args` | STRING |  | keyword arguments for the pca implementation, as key=value pairs in one quoted string, for example --pca_args="n_components=2". Values are read as Python literals where possible. |
| `--tsne_args` | STRING |  | keyword arguments for the tsne implementation, as key=value pairs in one quoted string, for example --tsne_args="n_components=2". Values are read as Python literals where possible. |
| `--umap_args` | STRING |  | keyword arguments for the umap implementation, as key=value pairs in one quoted string, for example --umap_args="n_components=2". Values are read as Python literals where possible. |
| `--scale` | NUMBER | 1000 | scale the largest axis so that the difference between the largest and smallest value is this. |
| `--three_d` |  |  | if set, then makes 3d projections. default: 2d projections. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
