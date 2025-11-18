"""Generate an enumerated report
Takes a genbank file which has been annotated using the domainator program, write tabulated descriptions of each contig, CDS, or domain.

write a tab-separated file
"""

import sys
from jsonargparse import ArgumentParser, ActionConfigFile
from domainator.utils import list_and_file_to_dict_keys, read_hmms
from domainator.select_by_cds import get_cds_neighborhood
from domainator import __version__, DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME, RawAndDefaultsFormatter
from pathlib import Path
from typing import List, Tuple, Union, Iterable, Dict
import json
from domainator.enum_report import DynamicArg, COLS_ARG_NAME, EnumTSVWriter, EnumHTMLWriter, get_analysis_names
import pyhmmer
import os
from io import IOBase


def read_hmms_to_iterators(hmm_files:Iterable[Union[str,os.PathLike, IOBase]]) -> Tuple[str, Iterable[pyhmmer.plan7.HMM]]:
    """
        hmm_files: a list of paths to .hmm files

        returns:
            a tuple of (name, iterator of pyhmmer.plan7.HMM)

    """
    out = list()
    for file in hmm_files:
        if isinstance(file, IOBase):
            name = file.name
        else:
            name = file
        name = os.path.basename(Path(name).stem)

        out.append((name, pyhmmer.plan7.HMMFile(file)))

    return out

def append_factory(spec):
    
    if len(spec) != 3:
        raise ValueError("append spec must be length 3")
    if spec[1] not in ["int", "float", "str"]:
        raise ValueError("append spec must be one of int, float, str")

    def append(source, rec):
        return spec[2]

    return {
        "columns": [spec[0]],
        "column_types": [spec[1]],
        "function": append,
    }


def hmmer_report(records:Tuple[str, Iterable[pyhmmer.plan7.Profile]], analyses, tsv_out_handle, html_out_handle, column_names, html_max_height):
    """
      input: 
        records: a list of tuples of tuple of (database_name, iterator of pyhmmer.plan7.HMM)
        by: One line in output for every by contig, cds or domain. default: by contig
        analyses: what data to report
        tsv_out_handle: a file handle to write the analysis table to.
        html_out_handle: a file handle to write html formatted output to.
        column_names: If supplied, then this list will be used instead of the default column names.
        html_max_height: Max height in pixels to set the html output to.
    """
    
    STATIC_ANALYSES = {"source": {"columns": ["source"], "column_types": ["str"], "function": lambda source,rec: source},
                "name": {"columns": ["name"], "column_types": ["str"], "function": lambda source,rec: "" if rec.name is None else rec.name.decode()},
                "acc": {"columns": ["acc"], "column_types": ["str"], "function": lambda source,rec: "" if rec.accession is None else rec.accession.decode()},
                "desc": {"columns": ["desc"], "column_types": ["str"], "function": lambda source,rec: "" if rec.description is None else rec.description.decode()},
                "length": {"columns": ["length"], "column_types": ["int"], "function": lambda source,rec: rec.M},
                "consensus": {"columns": ["consensus"], "column_types": ["str"], "function": lambda source,rec: "" if rec.consensus is None else rec.consensus},
                }
    DYNAMIC_ANALYSES = {"append": append_factory} # values are functions that return dicts of {"columns":[names_to_appear_in_output],  "column_types": [types_of_columns], "function": function taking rec, loc, tax as arguments and returning a scalar or list of scalars}
    
    headers = ["name"]
    column_types = ["str"]
    analyses_to_run = list()
    if analyses is None:
        analyses = list()
    for analysis in analyses:
        if isinstance(analysis, str): # if it's a string, then it's a static analysis
            headers.extend(STATIC_ANALYSES[analysis]["columns"])
            column_types.extend(STATIC_ANALYSES[analysis]["column_types"])
            analyses_to_run.append(STATIC_ANALYSES[analysis])
        else:
            generated_analysis = DYNAMIC_ANALYSES[analysis[0]](analysis[1])
            headers.extend(generated_analysis["columns"])
            column_types.extend(generated_analysis["column_types"])
            analyses_to_run.append(generated_analysis)

    
    if column_names is not None:
        if len(headers) != len(column_names):
            raise ValueError('The number of override column names must be the same as the number of columns.\nColumns:\n' + "; ".join(headers) + "\nOverride columns:\n" + "; ".join(column_names))
        headers = column_names

    writers = list()
    if tsv_out_handle is not None:
        writers.append(EnumTSVWriter(headers, column_types, tsv_out_handle))
    if html_out_handle is not None:
        writers.append(EnumHTMLWriter(headers, column_types, html_out_handle, html_max_height))

    for writer in writers:
        writer.write_header()
    

    for source, records in records:
        for rec in records:
            name = "" if rec.name is None else rec.name.decode()
            out_line = [name]
            for analysis in analyses_to_run:
                analysis_result = analysis["function"](source, rec)
                
                if isinstance(analysis_result, list) or isinstance(analysis_result, tuple): # doesn't handle generators
                    out_line.extend(analysis_result)
                else:
                    out_line.append(analysis_result)
            for writer in writers:
                writer.write_row(out_line)

    for writer in writers:
        writer.write_footer()

def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument('-i', '--input', nargs='+', required=False, default=None,
                        help="hmm filenames. If not supplied then reads from to stdin.")

    parser.add_argument('-o', '--output', default=None, required=False, 
                        help="The name of the output tsv file.")
    parser.add_argument('--html', default=None, required=False, 
                        help="Write a table in html format to this file.")
    parser.add_argument('--html_max_height', default=None, required=False, type=int,
                        help="Max height in pixels to set the html output to.")
   
    parser.add_argument('--column_names', nargs='+', default=None, required=False, type=str,
                        help="If supplied, then this list will be used instead of the default column names.")

    parser.add_argument('--source', action='append_const', dest=COLS_ARG_NAME, const="source",
                        help="report the name of the hmm file.")

    parser.add_argument('--acc', action='append_const', dest=COLS_ARG_NAME, const="acc",
                        help="report the start coordinate of the feature.")
    
    parser.add_argument('--desc', action='append_const', dest=COLS_ARG_NAME, const="desc",
                        help="return the end coordinate of the feature.")
    
    parser.add_argument('--length', action='append_const', dest=COLS_ARG_NAME, const="length",
                        help="return the length of the profile (i.e. the number of nodes).")
    
    parser.add_argument('--consensus', action='append_const', dest=COLS_ARG_NAME, const="consensus",
                        help="return the consensus residue line of the profile, if set..")
    
    parser.add_argument('--append', nargs=3, required=False, action=DynamicArg, dest=COLS_ARG_NAME, const="append",
                        help="Supply three strings, a column will be added with the first string as the column name, the second string as the column type (str, int, float) and the third string as the value for all rows.")
    
    parser.add_argument('--config', action=ActionConfigFile)

    params = parser.parse_args(argv)
   
    ### Figure out what input and output files ####
    if params.input is None:
        inputs = [sys.stdin]
    else:
        inputs = params.input
    
    out = None
    if params.output is not None:
        out = open(params.output, "w")

    html_out_handle = None
    if params.html is not None:
        html_out_handle = open(params.html, "w")

    if out is None and html_out_handle is None:
        raise ValueError("You must supply either --output or --html or both.")

    analyses = list()
    if hasattr(params, COLS_ARG_NAME): 
        analyses = getattr(params, COLS_ARG_NAME)

    # Run
    hmmer_report(read_hmms_to_iterators(inputs), analyses, out, html_out_handle, params.column_names, params.html_max_height)

    if params.output is not None:
        out.close()
    
    if params.html is not None:
        html_out_handle.close()

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
