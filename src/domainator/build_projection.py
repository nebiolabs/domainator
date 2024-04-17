"""generates 2d or 3d umap, tsne, or pca projections from dense or sparse matrices.

    output can be images, coordinates, or Cytoscape xgmml files.

"""
import argparse
import sys
from pathlib import Path
import warnings
import pandas as pd
import seaborn as sns
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA, TruncatedSVD
import umap
from matplotlib import pyplot as plt
import numpy as np
import json
from domainator import __version__, RawAndDefaultsFormatter
from domainator.data_matrix import DataMatrix
from domainator.cytoscape import write_cytoscape_xgmml
import scipy as sp
from domainator.color_genbank import read_color_table
from domainator.utils import get_palette


XCOLNAME="projection_X"
YCOLNAME="projection_Y"
ZCOLNAME="projection_Z"

#TODO: for projections that are less overweighted towards dense regions of sequence space, learn a PCA on the reduced sequence set then apply it to the expanded sequence set.




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
    
    data = data_matrix.data
    if data_matrix.sparse:
        data = data.toarray()

    return data




def build_projection(data_matrix: DataMatrix, algorithm, metadata_files, color_by, png_out, coords_out, xgmml, scale, three_dee=False, algorithm_parameters=None, color_table=None, color_table_out=None):
    """

        algorithm:
        metadata_files:
        color_by:
        png_out:
        coords_out:
        xgmml:
        scale:
        three_dee:
        algorithm_parameters: extra parameters to pass to the dimensionality reduction algorithm.
        color_table: dictionary mapping values in the color_by column to colors.
        color_table_out: tab separated file with two columns and no header, columns are: annotation, hex color. Written after the color table is updated with new colors, for example if using --color_by, but not supplying an external color table.

    """
    sparse = data_matrix.sparse
    n_components=2
    columns=[XCOLNAME,YCOLNAME]
    if three_dee:
        n_components=3
        columns += [ZCOLNAME]
    
    if algorithm_parameters is None:
        algorithm_parameters = dict()
    if algorithm == "tsne":
        if 'learning_rate' not in algorithm_parameters:
            algorithm_parameters["learning_rate"] = 'auto'
        if 'init' not in algorithm_parameters:
            algorithm_parameters["init"] = 'pca'
        data = data_matrix.data
        if sparse:
            if len(data_matrix.columns) > 50:
                #TODO: maybe warn that we're doing a dimenisionality reduction
                data = TruncatedSVD(n_components=50).fit_transform(data_matrix.data)
            else:
                data = data_matrix.data.toarray()
        if 'perplexity' not in algorithm_parameters:
            if data.shape[0] <= 30: # 30 is the default perplexity in scikit-learn tsne.
                algorithm_parameters["perplexity"] = data.shape[0] - 1
            else:
                algorithm_parameters["perplexity"] = 30 # TODO: what should the default be?
        embedded = TSNE(n_components=n_components, **algorithm_parameters).fit_transform(data)
    elif algorithm == "pca":
        if sparse:
            #TODO: maybe warn that we're doing TruncatedSVD instead of PCA
            embedded = TruncatedSVD(n_components=n_components).fit_transform(sp.sparse.csr_matrix(data_matrix.data)) #TODO: when sklearn is updated to support sparse arrays instead of matrices, delete the cast.
        else:
            embedded = PCA(n_components=n_components, **algorithm_parameters).fit_transform(data_matrix.data)
    elif algorithm == "umap":
        #doesn't matter if it's sparse or not umap will take it!
        if sparse:
            embedded = umap.UMAP(n_components=n_components, **algorithm_parameters).fit_transform(sp.sparse.csr_matrix(data_matrix.data)) #TODO: when sklearn is updated to support sparse arrays instead of matrices, delete the cast.
        else:    
            embedded = umap.UMAP(n_components=n_components, **algorithm_parameters).fit_transform(data_matrix.data)

    elif algorithm == "coords":
        embedded = check_coords(n_components, data_matrix,  **algorithm_parameters)

    embedded = center_and_scale(embedded, scale)

    projection = pd.DataFrame(embedded, columns=columns, index=data_matrix.rows)
    if coords_out is not None:
        projection.to_csv(coords_out,sep="\t")

    if metadata_files is not None:
        for file in metadata_files:
            metadata = pd.read_csv(file, sep="\t", index_col=0)
            projection = projection.merge(metadata,how="left", left_index=True, right_index=True)

    if color_by is not None:
        if color_table is None:
            color_table = get_palette(metadata[color_by].unique())

    if png_out is not None:
        if three_dee:
            warnings.warn("png output for 3d graphs not supported yet.")
        else:
            # change matplotlib backend
            old_backend = plt.get_backend()
            plt.switch_backend('agg') # this is necessary to run on a headless server
            sns.scatterplot(x=XCOLNAME,y=YCOLNAME, hue=color_by, data=projection, palette=color_table)
            plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
            plt.savefig(png_out, dpi=300, bbox_inches='tight')
            plt.clf()
            plt.switch_backend(old_backend)

    if xgmml is not None:
        write_cytoscape_xgmml(xgmml, projection, name=Path(xgmml).stem, x_col=XCOLNAME, y_col=YCOLNAME, z_col=ZCOLNAME, color_by=color_by, color_table=color_table)
    
    if color_table_out is not None:
        with open(color_table_out, "w") as out:
            for domain, color in color_table.items():
                out.write(f"{domain}\t{color}\n")


def main(argv):
    parser = argparse.ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-i', '--input', type=str, required=False,
                        help="The matrix that you want to project into 2D or 3D space, rows are observations (typically gene names), columns are features (typically similarity to other genes, or hmm profile scores). If passing in a coords matrix, use this (and not --sparse).")

    #TODO: add a way to train the projection on a subset of rows and/or columns ("training rows")
    #TODO: train just on reference sequences then project on all sequences
    #TODO: add a separate program for matrix analysis, such as density, non-zeros, histograms like they do in EFI
    #TODO: tree output (or as a separate script, or as part of the ssn script?), in the ssn 
    #TODO: auto-generate MSAs and hmms and logos and stuff from clusters

    parser.add_argument('--algorithm', type=str, required=False, default="pca", choices={"pca","tsne","umap", "coords"},
                        help="Which projection algorithm to use.")

    parser.add_argument('--metadata', type=str, nargs="+", required=False, default=None, 
                        help="tab separated files of sequence metadata.")
    
    parser.add_argument('--color_by', type=str, required=False, default=None,
                        help="Color the points in the output image based on this column of the metadata table.")
    parser.add_argument("--color_table", required=False, default=None, type=str, help="tab separated file with two columns and no header, columns are: annotation, hex color. For example: CCDB   cc0000")
    parser.add_argument("--color_table_out", required=False, default=None, type=str, help="tab separated file with two columns and no header, columns are: annotation, hex color. Written after the color table is updated with new colors, for example if using --color_by, but not supplying an external color table.")
    
    parser.add_argument('--png_out', type=str, required=False, default=None,
                        help="write a png file of the projection.")
    
    parser.add_argument('--xgmml', type=str, required=False, default=None,
                        help="write a Cytoscape xgmml file of the projection.")

    parser.add_argument('--coords_out', type=str, required=False, default=None,
                        help="write a tsv file with the projected coordinates.")

    # parser.add_argument('--gradient', type=str, required=False, default=None,
    #                     help="Color nodes based on a gradient instead of discrete value mapping.")

    parser.add_argument('--params', type=str, required=False, default=None,
                        help="additional parameters to pass to the algorithm, in the form of a json dict.")
    
    parser.add_argument('--scale', type=float, required=False, default=1000,
                        help="scale the largest axis so that the difference between the largest and smallest value is this. Default: 1000")
    
    parser.add_argument('--3d', action='store_true', default=False,
                        help="if set, then makes 3d projections. default: 2d projections.")
    
    #TODO: add a way to generate an interactive html output, could be with plotly. You could store the metadata as mouseover notes.


    params = parser.parse_args(argv)


    algorithm_parameters = {}
    if params.params:
        new_params = "{" + params.params + "}" 
        algorithm_parameters = json.loads(new_params)

    input_file = params.input

    color_table = None
    if params.color_table:
        color_table = read_color_table(params.color_table)

    # Run
    # def build_projection(data_matrix: DataMatrix, algorithm, metadata_files, color_by, png_out, coords_out, xgmml, scale, three_dee=False, algorithm_parameters=None)
    build_projection(data_matrix=DataMatrix.from_file(input_file), algorithm=params.algorithm, metadata_files=params.metadata, color_by=params.color_by, 
                     png_out=params.png_out, coords_out=params.coords_out, xgmml=params.xgmml, scale=params.scale, three_dee=vars(params)['3d'], 
                     algorithm_parameters=algorithm_parameters,
                     color_table=color_table, color_table_out=params.color_table_out)

def _entrypoint():
    main(sys.argv[1:])
    
if __name__ == '__main__':
    main(sys.argv[1:])
