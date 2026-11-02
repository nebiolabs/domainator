"""Build an HMM profile from a multiple sequence alignment (MSA)
    
    Allows the user to specify the ACC, NAME, and DESC fields of the HMM profile.
    Some of these options are not available in the standard hmmbuild tool from HMMER.
"""
from jsonargparse import ActionConfigFile
from pyhmmer.easel import MSAFile
from domainator.cli.common import add_input, add_output, check_outputs, declare_outputs, STDIN_PATHS, tool_parser
from domainator.hmmio import ALPHABET_NAMES, get_alphabet, is_nucleic_acid_alphabet, alphabet_name, open_writable_hmm_file
from domainator.cli.common import OUTPUT_COMPRESSION_HELP
import os
import shutil
import sys
import tempfile
import pyhmmer
from typing import Optional,BinaryIO,Union
import re

MSA_FORMATS_HELP = ("Acceptable formats are those of hmmbuild: Stockholm, Pfam, aligned FASTA, A2M, "
                    "Clustal, PHYLIP, SELEX and PSI-BLAST.")

def sanitize_string(s:str) -> str:
    return re.sub(r"[^ \w\d_\-\.;:]", "_", s)


def _spool(file:BinaryIO, directory:str) -> str:
    """Copy a stream into a file in `directory` and return the file's path.

    pyhmmer's MSAFile misbehaves on a file object: under pyhmmer 0.12 it segfaults the
    interpreter when the format cannot be guessed (empty, blank or unrecognized input),
    and a SELEX stream fails with "illegal residue(s)" once the alphabet has been
    guessed. Opened from a path it does neither. The copy has no suffix, so its format
    is guessed from the content alone.
    """
    path = os.path.join(directory, "msa")
    with open(path, "wb") as out:
        shutil.copyfileobj(file, out)
    return path


def _read_msa(path:str, alphabet:Optional[pyhmmer.easel.Alphabet]) -> pyhmmer.easel.MSA:
    """Read the first alignment of the MSA file at `path`.

    Opening by path lets Easel guess the format as hmmbuild does: from the first line,
    with the filename suffix (after a `.gz`) deciding what the content cannot (`.a2m`,
    `.phys`, `.pb`, `.slx`).
    """
    if os.path.getsize(path) == 0:
        raise ValueError("The input MSA is empty.")
    try:
        with MSAFile(path, digital=True, alphabet=alphabet) as msa_file:
            return msa_file.read()
    except ValueError as exc:
        if "Could not determine format" in str(exc):
            raise ValueError(f"Could not determine the format of the input MSA. {MSA_FORMATS_HELP}") from None
        if alphabet is None and "alphabet" in str(exc):
            raise ValueError(
                f"Could not determine the alphabet of the input MSA ({exc}). "
                f"Specify it explicitly with --alphabet ({', '.join(ALPHABET_NAMES)})."
            ) from None
        raise

def hmmer_build(file:Union[str,BinaryIO], *, alphabet:Optional[pyhmmer.easel.Alphabet]=None, name:Optional[str]=None, acc:Optional[str]=None, desc:Optional[str]=None, window_length:Optional[int]=None, window_beta:Optional[float]=None) -> pyhmmer.plan7.HMM:
    """Build a single plan7 profile from an MSA.

    Args:
        file: path to, or open binary handle on, an MSA in any format hmmbuild accepts.
            A path's suffix can decide the format (`.a2m`, `.phys`, `.pb`); a handle's
            format is guessed from its content.
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

    if isinstance(file, (str, os.PathLike)):
        msa = _read_msa(os.fspath(file), alphabet)
    else:
        with tempfile.TemporaryDirectory() as directory:
            msa = _read_msa(_spool(file, directory), alphabet)

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
    

def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, "Path of input msa. Acceptable formats are the same as for hmmbuild. The format is guessed from the content, and from the file suffix where the content is ambiguous (.a2m, .phys, .pb, .slx).", multiple=False)

    add_output(parser, "hmm", "Write the profile to this hmm file." + OUTPUT_COMPRESSION_HELP, primary=True)

    parser.add_argument("--name", default=None, required=True, type=str,
                            help="Name of the HMM profile.")
    parser.add_argument("--acc", default=None, required=False, type=str,
                            help="Accession of the HMM profile.")
    parser.add_argument("--desc", default=None, required=False, type=str,
                            help="Description of the HMM profile.")
    parser.add_argument("--alphabet", default=None, required=False, type=str.lower, choices=list(ALPHABET_NAMES),
                            help="Alphabet of the input MSA. If not supplied, inferred from the input, which can guess wrong on short or ambiguous alignments."
                        )
    parser.add_argument("--window_length", default=None, required=False, type=int,
                            help="Max expected hit length (the MAXL field of the profile), used by nhmmer for windowed E-values. Nucleotide alphabets only. If not supplied, hmmer picks one. Takes precedence over --window_beta, so only one of the two may be given.")
    parser.add_argument("--window_beta", default=None, required=False, type=float,
                            help="Tail mass from which to derive the window length. Nucleotide alphabets only. Only one of --window_length and --window_beta may be given.")

    parser.add_argument("--config", action=ActionConfigFile)
    declare_outputs(parser, ["hmm_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)


    check_outputs(parser, params)
    from_stdin = params.input in STDIN_PATHS
    input_file = sys.stdin.buffer if from_stdin else params.input
    output_handle = open_writable_hmm_file(params.hmm_out)

    alphabet = get_alphabet(params.alphabet)

    hmm = hmmer_build(file=input_file, alphabet=alphabet, name=params.name, acc=params.acc, desc=params.desc, window_length=params.window_length, window_beta=params.window_beta)
    hmm.write(output_handle)

    output_handle.close()
