"""Build an HMM profile from a multiple sequence alignment (MSA)
    
    Allows the user to specify the ACC, NAME, and DESC fields of the HMM profile.
    Some of these options are not available in the standard hmmbuild tool from HMMER.
"""
from jsonargparse import ArgumentParser, ActionConfigFile
from pyhmmer.easel import MSAFile
from domainator import __version__, RawAndDefaultsFormatter
from domainator.utils import ALPHABET_NAMES, get_alphabet, is_nucleic_acid_alphabet, alphabet_name
import os
import sys
import pyhmmer
from typing import Optional,BinaryIO,Union
import re

def sanitize_string(s:str) -> str:
    return re.sub(r"[^ \w\d_\-\.;:]", "_", s)


def _reject_empty_input(file) -> None:
    """Raise on empty input before pyhmmer can see it.

    MSAFile(io.BytesIO(b""), digital=True) segfaults the interpreter under
    pyhmmer 0.12.1, and that is reachable from the CLI through empty stdin.
    An empty path, by contrast, raises a normal error, so only streams need
    checking here.
    """
    if isinstance(file, (str, os.PathLike)):
        return
    peek = getattr(file, "peek", None)
    if peek is not None: # buffered readers, including sys.stdin.buffer
        if not peek(1):
            raise ValueError("The input MSA is empty.")
        return
    if getattr(file, "seekable", lambda: False)(): # BytesIO and friends
        position = file.tell()
        empty = not file.read(1)
        file.seek(position)
        if empty:
            raise ValueError("The input MSA is empty.")

def hmmer_build(file:Union[str,BinaryIO], alphabet:Optional[pyhmmer.easel.Alphabet]=None, name:Optional[str]=None, acc:Optional[str]=None, desc:Optional[str]=None, window_length:Optional[int]=None, window_beta:Optional[float]=None) -> pyhmmer.plan7.HMM:
    """Build a single plan7 profile from an MSA.

    Args:
        file: path to, or open binary handle on, an MSA in any format hmmbuild accepts.
        alphabet: alphabet of the MSA. If None, pyhmmer infers it from the residues,
            which can guess wrong on short or ambiguous alignments.
        name: value for the profile's NAME field.
        acc: value for the profile's ACC field.
        desc: value for the profile's DESC field.
        window_length: max expected hit length (the MAXL field), used by nhmmer for
            windowed E-values. Nucleotide alphabets only. Takes precedence over
            window_beta, so only one of the two may be given.
        window_beta: tail mass from which to derive the window length. Nucleotide
            alphabets only.

    Returns:
        the built HMM.
    """
    if window_length is not None and window_beta is not None:
        raise ValueError("Only one of window_length and window_beta may be specified, because window_length takes precedence over window_beta.")

    _reject_empty_input(file)

    try:
        with MSAFile(file, digital=True, alphabet=alphabet) as msa_file:
            msa = msa_file.read()
    except ValueError as exc:
        if alphabet is None and "alphabet" in str(exc):
            raise ValueError(
                f"Could not determine the alphabet of the input MSA ({exc}). "
                f"Specify it explicitly with --alphabet ({', '.join(ALPHABET_NAMES)})."
            ) from None
        raise

    if (window_length is not None or window_beta is not None) and not is_nucleic_acid_alphabet(msa.alphabet):
        raise ValueError(f"window_length and window_beta apply only to nucleotide alphabets, but the input MSA is {alphabet_name(msa.alphabet)}.")

    if name is not None:
        msa.name = sanitize_string(name).encode()
    elif msa.name is None:
        # Builder.build_msa() would otherwise fail with "Unable to name the HMM."
        raise ValueError("The input MSA is unnamed, so the profile cannot be named. Supply a name with --name.")
    if acc is not None:
        msa.accession = sanitize_string(acc).encode()
    if desc is not None:
        msa.description = sanitize_string(desc).encode()

    builder_kwargs = dict()
    if window_length is not None:
        builder_kwargs["window_length"] = window_length
    if window_beta is not None:
        builder_kwargs["window_beta"] = window_beta

    builder = pyhmmer.plan7.Builder(msa.alphabet, **builder_kwargs)
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
    parser.add_argument("--alphabet", default=None, required=False, type=str.lower, choices=set(ALPHABET_NAMES),
                            help="Alphabet of the input MSA. If not supplied, inferred from the input, which can guess wrong on short or ambiguous alignments."
                        )
    parser.add_argument("--window_length", default=None, required=False, type=int,
                            help="Max expected hit length (the MAXL field of the profile), used by nhmmer for windowed E-values. Nucleotide alphabets only. If not supplied, hmmer picks one. Takes precedence over --window_beta, so only one of the two may be given.")
    parser.add_argument("--window_beta", default=None, required=False, type=float,
                            help="Tail mass from which to derive the window length. Nucleotide alphabets only. Only one of --window_length and --window_beta may be given.")

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

    alphabet = get_alphabet(params.alphabet)

    hmm = hmmer_build(file=input_file, alphabet=alphabet, name=params.name, acc=params.acc, desc=params.desc, window_length=params.window_length, window_beta=params.window_beta)
    hmm.write(output_handle)

    if params.input is not None:
        input_file.close()

    if params.output is not None:
        output_handle.close()

def _entrypoint():
    main(sys.argv[1:])

if __name__ == "__main__":
    main(sys.argv[1:])