"""
    Script to build a HMM profile from a multiple sequence alignment (MSA).
    Allows the user to specify the ACC, NAME, and DESC fields of the HMM profile.    
"""
from jsonargparse import ArgumentParser, ActionConfigFile
from pyhmmer.easel import MSAFile
from domainator import __version__, RawAndDefaultsFormatter
import sys
import pyhmmer
from typing import Optional,BinaryIO,Union
import re

def sanitize_string(s:str) -> str:
    return re.sub("[^ \w\d_\-\.;:]", "_", s)

def hmmer_build(file:Union[str,BinaryIO], alphabet:Optional[pyhmmer.easel.Alphabet]=None, name:Optional[str]=None, acc:Optional[str]=None, desc:Optional[str]=None) -> pyhmmer.plan7.HMM:
    with MSAFile(file, digital=True, alphabet=alphabet) as msa_file:
        msa = msa_file.read()

    msa.name = sanitize_string(name).encode()
    if acc is not None:
        msa.accession = sanitize_string(acc).encode()
    if desc is not None:
        msa.description = sanitize_string(desc).encode()
    
    builder = pyhmmer.plan7.Builder(msa.alphabet)
    background = pyhmmer.plan7.Background(msa.alphabet)
    hmm, _, _ = builder.build_msa(msa, background)
    return hmm
    

def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument("-i", "--input", default=None, required=False, type=str,
                          help="Path of input msa. If not supplied, reads from stdin. Acceptable formats are the same as for hmmbuild.")

    parser.add_argument("-o", "--output", default=None, required=False,  type=str,
                        help="hmm output file path. If not supplied writes to stdout.")

    parser.add_argument("--name", default=None, required=True, type=str,
                            help="Name of the HMM profile.")
    parser.add_argument("--acc", default=None, required=False, type=str,
                            help="Accession of the HMM profile.")
    parser.add_argument("--desc", default=None, required=False, type=str,
                            help="Description of the HMM profile.")
    parser.add_argument("--alphabet", default=None, required=False, type=str.lower, choices={"amino", "dna", "rna"},)

    parser.add_argument("--config", action=ActionConfigFile)

    params = parser.parse_args(argv)


    if params.input is None:
        input_file = sys.stdin.buffer
    else:
        input_file = open(params.input, "rb")

    if params.output is None:
        output_handle = sys.stdout.buffer
    else:
        output_handle = open(params.output, "wb")

    alphabet = None
    if params.alphabet == "amino":
        alphabet = pyhmmer.easel.Alphabet.amino()
    elif params.alphabet == "dna":
        alphabet = pyhmmer.easel.Alphabet.dna()
    elif params.alphabet == "rna":
        alphabet = pyhmmer.easel.Alphabet.rna()

    hmm = hmmer_build(file=input_file, alphabet=alphabet, name=params.name, acc=params.acc, desc=params.desc)
    hmm.write(output_handle)

    if params.input is not None:
        input_file.close()

    if params.output is not None:
        output_handle.close()

def _entrypoint():
    main(sys.argv[1:])

if __name__ == "__main__":
    main(sys.argv[1:])