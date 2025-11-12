""" reads a DataMatrix and generates a report, including a histogram of the non-zero values in the matrix.

"""

#TODO: deduplicate symmetric matrices. For example, only take the values from the lower triangualr matrix and ignore the upper triangular matrix of symmetric matrices.
#TODO: ignore the diagonal of symmetric matrices.
#TODO: add cumulative distribution plot for example to help picking a threshold with a convenient number of edges, or clusters

import argparse
import sys

from domainator.data_matrix import DataMatrix
from domainator import __version__, RawAndDefaultsFormatter
from bashplotlib.histogram import plot_hist
from bashplotlib.scatterplot import plot_scatter
import io
from contextlib import redirect_stdout
from domainator.summary_report import histplot_base64
import statistics

class SummaryTextWriter():
    def __init__(self, out_handle): 
        self.out_handle = out_handle
    
    def write_header(self, matrix, values):
        print(f"""Matrix Report
Total values: {matrix.size}
Non-zero values: {len(values)}
Mean: {statistics.mean(values)}
Median: {statistics.median(values)}
Min: {min(values)}
Max: {max(values)}
        """, file=self.out_handle)
        
    def write_hist(self, values):
        if len(values) != 0:
            hist_str = io.StringIO()
            with redirect_stdout(hist_str):
                plot_hist(values, bincount=50, title="Matrix Values", height=20.0, xlab=True, regular=True, showSummary=False)
            hist_str = hist_str.getvalue() # get the string
            hist_str = hist_str.replace('[39m', '') # remove color codes
        else:
            hist_str = ""
        print(f"""
{hist_str}""", file=self.out_handle)
    def write_footer(self):
        pass

class SummaryHTMLWriter():
    def __init__(self, out_handle): 
        self.out_handle = out_handle
    
    def write_header(self, matrix, values):
        print(f"""<!doctype html><html><head><meta charset="UTF-8" /><title>Matrix Report</title></head><body>""", 
              file=self.out_handle)
        print(f"""<h1>Matrix Report</h1>""", file=self.out_handle)
        print(f"""<div><h2>Summary</h2>""", file=self.out_handle)
        print(f"""<table>""", file=self.out_handle)
        print(f"""<tr><th>Total values</th><td>{matrix.size}</td></tr>""", file=self.out_handle)
        print(f"""<tr><th>Non-zero values</th><td>{len(values)}</td></tr>""", file=self.out_handle)
        print(f"""<tr><th>Mean</th><td>{statistics.mean(values) : .1f}</td></tr>""", file=self.out_handle)
        print(f"""<tr><th>Median</th><td>{statistics.median(values) : .1f}</td></tr>""", file=self.out_handle)
        print(f"""<tr><th>Min</th><td>{min(values) : .1f}</td></tr>""", file=self.out_handle)
        print(f"""<tr><th>Max</th><td>{max(values) : .1f}</td></tr>""", file=self.out_handle)
        print(f"""</table>""", file=self.out_handle)
        print(f"""</div>""", file=self.out_handle)
    def write_hist(self, values):
        print(f"""<div><h2>Matrix Values</h2>""", file=self.out_handle)
        print(f"""<img alt="Contig lengths histogram" src="data:image/png;base64, {histplot_base64(values)}" />""", file=self.out_handle)
        print(f"""</div>""", file=self.out_handle)

    def write_footer(self):
        print("</body></html>", file=self.out_handle)


def matrix_report(matrix:DataMatrix, out_text_handle, out_html_handle):
    """
        Write a report on the matrix to the given handles.
    """
    longform_outputs = list()
    if out_text_handle is not None:
        longform_outputs.append(SummaryTextWriter(out_text_handle))
    if out_html_handle is not None:
        longform_outputs.append(SummaryHTMLWriter(out_html_handle))
    
    values = list()
    for (row_name, column_name, value) in matrix.iter_data():
        if value != 0:
            values.append(value)
    values.sort()
    for output in longform_outputs:
        output.write_header(matrix, values)
        output.write_hist(values)
        output.write_footer()

def main(argv):
    parser = argparse.ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument("-i", "--input", default=None, required=True, type=str, help="name of input matrix file.")
    
    parser.add_argument('-o', '--output', default=None, required=False,
                        help="Text file to write output to.")
    
    parser.add_argument('--html', default=None, required=False,
                        help="html file to write output to.")
    
    params = parser.parse_args(argv)

    if params.output is None and params.html is None:
        parser.error("No output file specified. Use --output and/or --html to specify at least one output file.")


    out = None
    if params.output is not None:
        out = open(params.output, "w")
    
    out_html_handle = None
    if params.html is not None:
        out_html_handle = open(params.html, "w")
 

    ### Run
    matrix_report(
        DataMatrix.from_file(params.input),
        out,
        out_html_handle,
    )

    if params.output is not None:
        out.close()
    
    if params.html is not None:
        out_html_handle.close()

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':

    main(sys.argv[1:])
  
