"""Writes record-wise information from hmm files into tab-separated or html tables

Amino acid, DNA, and RNA profiles are all supported; use --alphabet to report which
alphabet each profile uses.
"""

import contextlib
from jsonargparse import ActionConfigFile
from domainator.cli.common import add_input, add_output, check_outputs, declare_outputs, tool_parser
from domainator.outputs import output_file
from domainator.hmmio import alphabet_name, iter_hmms_with_alphabet
from domainator.seqio import db_name_from_path
from typing import Tuple, Union, Iterable
from domainator.report_columns import Column, add_columns_argument, parse_columns
from domainator.table_writers import HTMLTableWriter, JSONTableWriter, TSVTableWriter
import pyhmmer
import os
from io import IOBase



COLUMNS = (
    Column("source", "the name of the hmm file."),
    Column("acc", "the profile's accession."),
    Column("desc", "the profile's description."),
    Column("length", "the length of the profile (the number of nodes)."),
    Column("consensus", "the profile's consensus residue line, if set."),
    Column("alphabet", "the profile's alphabet (amino, DNA, or RNA)."),
    Column("append", "a column with this name and type (str, int or float) holding this value in every row.", ("name", "type", "value")),
)

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
        name = db_name_from_path(name)

        _alphabet, profiles = iter_hmms_with_alphabet(file, role="input")
        out.append((name, profiles))

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


def hmmer_report(records:Tuple[str, Iterable[pyhmmer.plan7.Profile]], *, analyses, tsv_out_handle, html_out_handle, column_names, html_max_height, json_out_handle=None):
    """
      input:
        records: a list of tuples of tuple of (database_name, iterator of pyhmmer.plan7.HMM)
        by: One line in output for every by contig, cds or domain. default: by contig
        analyses: what data to report
        tsv_out_handle: a file handle to write the analysis table to.
        html_out_handle: a file handle to write html formatted output to.
        column_names: If supplied, then this list will be used instead of the default column names.
        html_max_height: Max height in pixels to set the html output to.
        json_out_handle: a file handle to write streaming NDJSON (one object per profile) to.
    """
    
    STATIC_ANALYSES = {"source": {"columns": ["source"], "column_types": ["str"], "function": lambda source,rec: source},
                "name": {"columns": ["name"], "column_types": ["str"], "function": lambda source,rec: "" if rec.name is None else rec.name},
                "acc": {"columns": ["acc"], "column_types": ["str"], "function": lambda source,rec: "" if rec.accession is None else rec.accession},
                "desc": {"columns": ["desc"], "column_types": ["str"], "function": lambda source,rec: "" if rec.description is None else rec.description},
                "length": {"columns": ["length"], "column_types": ["int"], "function": lambda source,rec: rec.M},
                "consensus": {"columns": ["consensus"], "column_types": ["str"], "function": lambda source,rec: "" if rec.consensus is None else rec.consensus},
                "alphabet": {"columns": ["alphabet"], "column_types": ["str"], "function": lambda source,rec: alphabet_name(rec.alphabet)},
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
        writers.append(TSVTableWriter(headers, column_types, tsv_out_handle))
    if html_out_handle is not None:
        writers.append(HTMLTableWriter(headers, column_types, html_out_handle, html_max_height))
    if json_out_handle is not None:
        writers.append(JSONTableWriter(headers, column_types, json_out_handle))

    for writer in writers:
        writer.write_header()
    

    for source, records in records:
        for rec in records:
            name = "" if rec.name is None else rec.name
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

def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, "hmm filenames.")

    add_output(parser, "tsv", "Write the table to this tab-separated file.", primary=True)
    add_output(parser, "html", "Write the table in html format to this file.")
    add_output(parser, "json", "Write streaming NDJSON (one compact JSON object per profile) to this file.")
    parser.add_argument('--html_max_height', default=None, required=False, type=int,
                        help="Max height in pixels to set the html output to.")

    parser.add_argument('--column_names', nargs='+', default=None, required=False, type=str,
                        help="If supplied, then this list will be used instead of the default column names.")

    add_columns_argument(parser, COLUMNS)

    
    
    

    
    
    parser.add_argument('--config', action=ActionConfigFile)
    declare_outputs(parser, ["tsv_out", "html_out", "json_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
   
    check_outputs(parser, params)

    analyses = parse_columns(params.columns or [], COLUMNS)

    # Run
    with contextlib.ExitStack() as outputs:
        def open_output(path):
            return None if path is None else outputs.enter_context(output_file(path))
        hmmer_report(read_hmms_to_iterators(params.input), analyses=analyses, tsv_out_handle=open_output(params.tsv_out), html_out_handle=open_output(params.html_out),
                     column_names=params.column_names, html_max_height=params.html_max_height, json_out_handle=open_output(params.json_out))
