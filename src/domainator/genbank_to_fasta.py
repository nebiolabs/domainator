"""
    Convert a GenBank file to a FASTA file.
"""
from jsonargparse import ArgumentParser, ActionConfigFile
import sys
from domainator.Bio import SeqIO
from domainator import __version__, RawAndDefaultsFormatter
from domainator.utils import parse_seqfiles

# TODO: replace this with a more general purpose script that can convert between any of the supported formats
# allow specification of various paramters, such as species, accession, etc.

def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument('-i', '--input', nargs='+', required=False, default=None,
                       help="Genbank filenames. If not supplied, reads from stdin.")

    parser.add_argument('-o', '--output', default=None, required=False,
                        help="the name of the output fasta file. If not supplied writes to stdout.")

    parser.add_argument('--config', action=ActionConfigFile)

    params = parser.parse_args(argv)

    ### Figure out what input and output files ####
    
    if params.input is None:
        genbanks = [sys.stdin]
    else:
        genbanks = params.input

    if params.output is None:
        out = sys.stdout
    else:
        out = open(params.output, "w")

    for rec in parse_seqfiles(genbanks, filetype_override="genbank"):
        SeqIO.write(rec, out, "fasta")
    
    if params.output is not None:
        out.close()

def _entrypoint():
    main(sys.argv[1:])

if __name__ == "__main__":
    main(sys.argv[1:])
