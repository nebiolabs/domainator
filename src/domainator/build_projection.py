"""Generates 2d or 3d umap, tsne, or pca projections from matrices

    output can be an interactive HTML plot, coordinates, or Cytoscape xgmml files.

"""
from typing import Any, Dict, List, Optional
import inspect
from jsonargparse import ActionConfigFile
from pathlib import Path
import pandas as pd
import numpy as np
from domainator.cli.common import add_input, add_output, check_outputs, declare_outputs, parse_keyword_args, tool_parser
from domainator.data_matrix import DataMatrix
from domainator.external import import_optional
from domainator.cytoscape import write_cytoscape_xgmml
import scipy as sp
from domainator.colors import get_palette, read_color_table
from domainator.viewer_bundle import read_node_metadata, viewer_bundle, viewer_defaults
from domainator.viewer_html import write_viewer_outputs


XCOLNAME="projection_X"
YCOLNAME="projection_Y"
ZCOLNAME="projection_Z"

DENSE_PROJECTION_MIN_FILL = 0.5
"""Matrices with at least this fraction of non-zero cells are projected as dense arrays."""

#TODO: for projections that are less overweighted toward dense regions of sequence space, learn a PCA on the reduced sequence set then apply it to the expanded sequence set.




def center_and_scale(matrix, max_value):
    """
        matrix should be dimensions, pts, projection_dimension
        moves the smallest values along each axis to zero (so that all coordinates are positive)
        Scales the axes proportionally so that the largest value on either axis is max_value
    """

    mins = np.min(matrix, axis=0)
    for i in range(matrix.shape[1]):
        matrix[:,i] -= mins[i]
    matrix = matrix * ( max_value / matrix.max())
    return matrix

def check_coords(n_components, data_matrix: DataMatrix):
    mat_columns = data_matrix.columns
    if n_components == 2:
        if len(mat_columns) != 2:
            raise ValueError(f"Input matrix should have two columns, not: {mat_columns}")
        if not ((mat_columns[0] == XCOLNAME) and (mat_columns[1] == YCOLNAME)):
            raise ValueError(f"Input matrix should have columns should be {[XCOLNAME,YCOLNAME]}, not: {mat_columns}")
    
    if n_components == 3:
        if len(mat_columns) != 3:
            raise ValueError(f"Input matrix should have three columns, not: {mat_columns}")
        if not ((mat_columns[0] == XCOLNAME) and (mat_columns[1] == YCOLNAME) and (mat_columns[2] == ZCOLNAME)):
            raise ValueError(f"Input matrix should have columns should be (in order) {[XCOLNAME,YCOLNAME,ZCOLNAME]}, not: {mat_columns}")
    
    return data_matrix.toarray()


def build_projection(data_matrix: DataMatrix, *, algorithm: str = "pca", metadata: Optional[List[str]] = None,
                     color_by: Optional[str] = None, html_out: Optional[str] = None, dvb_out: Optional[str] = None,
                     coords_out: Optional[str] = None,
                     xgmml_out: Optional[str] = None, scale: Optional[float] = 1000,
                     three_d: bool = False, algorithm_parameters: Optional[Dict[str, Any]] = None,
                     color_table: Optional[Dict[str, str]] = None, color_table_out: Optional[str] = None,
                     palette_offset: int = 0, label_by: Optional[str] = None, categorical: Optional[List[str]] = None,
                     name: Optional[str] = None) -> None:
    """Project the rows of ``data_matrix`` into 2D or 3D and write the projection.

    Args:
        data_matrix: rows are the points (typically sequences), columns their features
            (typically scores against other sequences or profiles).
        algorithm: "pca", "tsne", "umap", or "coords" (the matrix already holds the
            coordinates).
        metadata: tab-separated metadata files, joined to the points by their first column.
        color_by: color the points by this metadata column.
        html_out: write the Domainator Viewer, showing the projection, to this HTML file.
        dvb_out: write the viewer bundle (gzip-compressed JSON, .dvb) to this file.
        coords_out: write the coordinates to this tab-separated file.
        xgmml_out: write the projection to this Cytoscape XGMML file.
        scale: scale the coordinates so that the widest axis spans this much.
        three_d: project into 3D instead of 2D.
        algorithm_parameters: keyword arguments for the projection algorithm. For "pca" on a sparse
            matrix, which is projected with scikit-learn's TruncatedSVD instead of PCA, only the
            keywords the two share apply.
        color_table: maps ``color_by`` values to colors. None: assign colors automatically.
        color_table_out: write the color table, with any colors assigned, to this path: two
            tab-separated columns (annotation, hex color) and no header.
        palette_offset: skip this many colors at the start of the palette when assigning
            colors to color_by values. Ignored when color_table is given.
        label_by: the metadata column the viewer labels points by at first.
        categorical: numeric metadata columns the viewer colors as categories, not a gradient.
        name: the projection's name, shown in the viewer. Defaults to the html or dvb output's stem.
    """
    if color_table_out is not None and color_by is None and color_table is None:
        raise ValueError("--color_table_out requires either --color_by or --color_table, otherwise there are no colors to write.")
    # Mostly-filled matrices (e.g. profile scores, distances) are projected densely; mostly-empty
    # ones (e.g. similarity networks) stay sparse, which changes PCA to TruncatedSVD.
    is_sparse = data_matrix.nnz < DENSE_PROJECTION_MIN_FILL * data_matrix.size
    n_components=2
    columns=[XCOLNAME,YCOLNAME]
    if three_d:
        n_components=3
        columns += [ZCOLNAME]
    
    if algorithm_parameters is None:
        algorithm_parameters = dict()
    if algorithm == "tsne":
        from sklearn.decomposition import TruncatedSVD
        from sklearn.manifold import TSNE
        if 'learning_rate' not in algorithm_parameters:
            algorithm_parameters["learning_rate"] = 'auto'
        if 'init' not in algorithm_parameters:
            algorithm_parameters["init"] = 'pca'
        if is_sparse and len(data_matrix.columns) > 50:
            #TODO: maybe warn that we're doing a dimenisionality reduction
            data = TruncatedSVD(n_components=50).fit_transform(data_matrix.data)
        else:
            data = data_matrix.toarray()
        if 'perplexity' not in algorithm_parameters:
            if data.shape[0] <= 30: # 30 is the default perplexity in scikit-learn tsne.
                algorithm_parameters["perplexity"] = data.shape[0] - 1
            else:
                algorithm_parameters["perplexity"] = 30 # TODO: what should the default be?
        embedded = TSNE(n_components=n_components, **algorithm_parameters).fit_transform(data)
    elif algorithm == "pca":
        from sklearn.decomposition import PCA, TruncatedSVD
        if is_sparse:
            #TODO: maybe warn that we're doing TruncatedSVD instead of PCA
            unsupported = sorted(set(algorithm_parameters) - set(inspect.signature(TruncatedSVD).parameters))
            if unsupported:
                raise ValueError(f"A sparse matrix is projected with TruncatedSVD instead of PCA, which does not take {', '.join(unsupported)}.")
            embedded = TruncatedSVD(n_components=n_components, **algorithm_parameters).fit_transform(sp.sparse.csr_matrix(data_matrix.data)) #TODO: when sklearn is updated to support sparse arrays instead of matrices, delete the cast.
        else:
            embedded = PCA(n_components=n_components, **algorithm_parameters).fit_transform(data_matrix.toarray())
    elif algorithm == "umap":
        umap = import_optional("umap", 'Install it with: pip install "domainator[projection]" (or pip install umap-learn).')
        #doesn't matter if it's sparse or not umap will take it!
        if is_sparse:
            embedded = umap.UMAP(n_components=n_components, **algorithm_parameters).fit_transform(sp.sparse.csr_matrix(data_matrix.data)) #TODO: when sklearn is updated to support sparse arrays instead of matrices, delete the cast.
        else:    
            embedded = umap.UMAP(n_components=n_components, **algorithm_parameters).fit_transform(data_matrix.toarray())

    elif algorithm == "coords":
        if algorithm_parameters:
            raise ValueError(f"algorithm parameters do not apply to --algorithm coords, which reads the coordinates as given: {algorithm_parameters}")
        embedded = check_coords(n_components, data_matrix)

    embedded = center_and_scale(embedded, scale)

    projection = pd.DataFrame(embedded, columns=columns, index=data_matrix.rows)
    if coords_out is not None:
        projection.to_csv(coords_out,sep="\t")

    node_data = read_node_metadata(metadata, data_matrix.rows)
    projection = projection.merge(node_data, how="left", left_index=True, right_index=True)

    given_color_table = color_table
    if color_by is not None:
        if color_table is None:
            color_table = get_palette(projection[color_by].unique(), offset=palette_offset)

    if html_out is not None or dvb_out is not None:
        if name is None:
            output = html_out if html_out is not None else dvb_out
            name = "projection" if hasattr(output, "write") else Path(output).stem
        # The viewer assigns default colors as get_palette does; an offset palette goes in as a table.
        viewer_colors = given_color_table if given_color_table is not None or palette_offset == 0 else color_table
        defaults = viewer_defaults(node_data, color_by=color_by, color_table=viewer_colors, label_by=label_by, categorical=categorical)
        view = {"algorithm": algorithm, "axes": columns, "coords": [projection[axis].tolist() for axis in columns]}
        write_viewer_outputs(viewer_bundle(name, data_matrix.rows, node_data, defaults, {"projection": view}),
                             html_out=html_out, dvb_out=dvb_out)

    if xgmml_out is not None:
        write_cytoscape_xgmml(xgmml_out, projection, name=Path(xgmml_out).stem, x_col=XCOLNAME, y_col=YCOLNAME, z_col=ZCOLNAME, color_by=color_by, color_table=color_table)
    
    if color_table_out is not None:
        with open(color_table_out, "w") as out:
            for domain, color in color_table.items():
                out.write(f"{domain}\t{color}\n")


def build_parser():
    parser = tool_parser(__doc__)
    
    add_input(parser, "The matrix that you want to project into 2D or 3D space, rows are observations (typically gene names), columns are features (typically similarity to other genes, or hmm profile scores). Format can be tab-separated text or Domainator hdf5.", multiple=False, stdin=False)

    #TODO: add a way to train the projection on a subset of rows and/or columns ("training rows")
    #TODO: train just on reference sequences then project on all sequences
    #TODO: add a separate program for matrix analysis, such as density, non-zeros, histograms like they do in EFI
    #TODO: tree output (or as a separate script, or as part of the ssn script?), in the ssn 
    #TODO: auto-generate MSAs and hmms and logos and stuff from clusters

    parser.add_argument('--algorithm', type=str, required=False, default="pca", choices=["pca", "tsne", "umap", "coords"],
                        help="Which projection algorithm to use.")

    parser.add_argument('--metadata', type=str, nargs="+", required=False, default=None, 
                        help="tab separated files of sequence metadata.")
    
    parser.add_argument('--color_by', type=str, required=False, default=None,
                        help="Color the points in the HTML plot and the xgmml by this column of the metadata table.")
    parser.add_argument("--color_table", required=False, default=None, type=str, help="tab separated file with two columns and no header, columns are: annotation, hex color. For example: CCDB   cc0000")
    parser.add_argument("--color_table_out", required=False, default=None, type=str, help="tab separated file with two columns and no header, columns are: annotation, hex color. Written after the color table is updated with new colors, for example if using --color_by, but not supplying an external color table.")
    parser.add_argument("--palette_offset", required=False, default=0, type=int, help="Skip this many colors at the start of the palette before assigning colors to --color_by values. For example, if a plot colored by one column used the first 5 colors, --palette_offset 5 colors another column without reusing them. Ignored with --color_table.")
    parser.add_argument("--label_by", type=str, default=None,
                        help="Label the points in the html viewer by this metadata column at first.")
    parser.add_argument("--categorical", type=str, nargs="+", required=False, default=None,
                        help="Numeric metadata columns (for example integer cluster numbers) that the html viewer should color as discrete categories instead of a gradient.")
    parser.add_argument("--name", type=str, default=None,
                        help="The projection's name, shown in the html viewer. Defaults to the html (or dvb) output's file name.")
    
    add_output(parser, "xgmml", "Write the projection to this Cytoscape xgmml file.", primary=True, stdout=False)
    add_output(parser, "html", "Write the Domainator Viewer, showing the projection, to this HTML file. It selects points into a metadata table and saves pictures as SVG or PNG.", stdout=False)
    add_output(parser, "dvb", "Write the gzip-compressed JSON Domainator Viewer bundle (.dvb) to this file. The viewer opens it, and saves sessions in the same format.", stdout=False)
    add_output(parser, "coords", "Write the projected coordinates to this tab-separated file.", stdout=False)

    for algorithm in ("pca", "tsne", "umap"):
        parser.add_argument(f'--{algorithm}_args', type=str, required=False, default=None,
                            help=f"keyword arguments for the {algorithm} implementation, as key=value pairs in one quoted string, "
                                 f"for example --{algorithm}_args=\"n_components=2\". Values are read as Python literals where possible.")
    
    parser.add_argument('--scale', type=float, required=False, default=1000,
                        help="scale the largest axis so that the difference between the largest and smallest value is this.")
    
    parser.add_argument('--three_d', action='store_true', default=False,
                        help="if set, then makes 3d projections. default: 2d projections.")
    
    parser.add_argument('--config', action=ActionConfigFile)

    declare_outputs(parser, ["xgmml_out", "html_out", "dvb_out", "coords_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
    check_outputs(parser, params)

    for algorithm in ("pca", "tsne", "umap"):
        if algorithm != params.algorithm and getattr(params, f"{algorithm}_args"):
            parser.error(f"--{algorithm}_args was given, but --algorithm is {params.algorithm}.")
    algorithm_args = getattr(params, f"{params.algorithm}_args", None)
    algorithm_parameters = parse_keyword_args(algorithm_args, f"--{params.algorithm}_args")

    input_file = params.input

    color_table = None
    if params.color_table:
        color_table = read_color_table(params.color_table)

    # Run
    build_projection(data_matrix=DataMatrix.from_file(input_file), algorithm=params.algorithm, metadata=params.metadata, color_by=params.color_by, 
                     html_out=params.html_out, dvb_out=params.dvb_out, coords_out=params.coords_out, xgmml_out=params.xgmml_out, scale=params.scale, three_d=params.three_d, 
                     algorithm_parameters=algorithm_parameters,
                     color_table=color_table, color_table_out=params.color_table_out,
                     palette_offset=params.palette_offset, label_by=params.label_by, categorical=params.categorical,
                     name=params.name)
