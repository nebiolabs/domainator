"""
    Convert a GenBank file to a FASTA file
"""
from typing import Iterable

from jsonargparse import ActionConfigFile

from domainator.cli.common import add_input, add_output, check_outputs, declare_outputs, tool_parser
from domainator.outputs import write_seq_outputs
from domainator.seqio import parse_seqfiles


def genbank_to_fasta(input: Iterable[str], *, fasta_out: str) -> None:
    """Write the records of the GenBank files ``input`` to the FASTA file ``fasta_out`` (``-`` for standard output)."""
    write_seq_outputs(parse_seqfiles(input, filetype_override="genbank"), fasta=fasta_out)


def build_parser():
    parser = tool_parser(__doc__)
    add_input(parser, 'Genbank filenames.')
    add_output(parser, "fasta", "Write the records to this FASTA file.", primary=True)
    parser.add_argument('--config', action=ActionConfigFile)
    declare_outputs(parser, ["fasta_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
    check_outputs(parser, params)
    genbank_to_fasta(params.input, fasta_out=params.fasta_out)
