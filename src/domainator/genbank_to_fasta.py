"""
    Convert a GenBank file to a FASTA file
"""
from jsonargparse import ActionConfigFile

from domainator.cli.common import add_input, add_output, check_outputs, tool_parser
from domainator.outputs import write_seq_outputs
from domainator.seqio import parse_seqfiles


def build_parser():
    parser = tool_parser(__doc__)
    add_input(parser, 'Genbank filenames.')
    add_output(parser, "fasta", "Write the records to this FASTA file.", primary=True)
    parser.add_argument('--config', action=ActionConfigFile)
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
    check_outputs(parser, params, ["fasta_out"])
    write_seq_outputs(parse_seqfiles(params.input, filetype_override="genbank"), fasta=params.fasta_out)
