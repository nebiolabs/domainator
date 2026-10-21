"""Reading and writing HMMER profiles, Infernal models and reference sequence sets.
"""
from collections import OrderedDict
from io import IOBase
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple, Union
import gzip
import math
import os
import sys
import warnings

import pyhmmer
import pyinfernal.cm

from domainator.seqio import STDIN_PATHS, STDOUT_PATHS, _UnclosedStdout, db_name_from_path, detect_compression, is_compressed_path


ALPHABET_NAMES = ("amino", "dna", "rna")
"""Accepted values for the ``--alphabet`` arguments of the hmmer_* tools."""

_ALPHABET_FACTORIES = {
    "amino": pyhmmer.easel.Alphabet.amino,
    "dna": pyhmmer.easel.Alphabet.dna,
    "rna": pyhmmer.easel.Alphabet.rna,
}

_AMINO_ALPHABET_K = 20
"""Size of the amino acid alphabet. The profile-profile scoring thresholds in
domainator hmmer_search were calibrated against it, so it is the reference for
alphabet_score_scale()."""


def get_alphabet(name: Optional[str]) -> Optional[pyhmmer.easel.Alphabet]:
    """Convert an --alphabet argument into a pyhmmer Alphabet.

    Args:
        name: one of "amino", "dna", "rna" (case insensitive), or None.

    Returns:
        The corresponding Alphabet, or None if name is None, which callers
        interpret as "infer the alphabet from the input".
    """
    if name is None:
        return None
    try:
        return _ALPHABET_FACTORIES[name.lower()]()
    except KeyError:
        raise ValueError(
            f"Unknown alphabet: {name}. Must be one of: {', '.join(ALPHABET_NAMES)}."
        ) from None


def is_nucleic_acid_alphabet(alphabet: pyhmmer.easel.Alphabet) -> bool:
    """True for DNA and RNA alphabets, False for amino."""
    return alphabet.is_dna() or alphabet.is_rna()


def alphabet_name(alphabet: pyhmmer.easel.Alphabet) -> str:
    """Short display name of an alphabet: "amino", "DNA", or "RNA"."""
    return alphabet.type


def alphabet_score_scale(alphabet: pyhmmer.easel.Alphabet) -> float:
    """Factor for rescaling protein-calibrated profile-profile scores to another alphabet.

    The per-column score computed by hmmer_search.Saa is
    log(sum_i(q_i * r_i / background_i)), which for a perfectly conserved column
    over a uniform background tops out at log(K). The score cutoffs and alignment
    midline thresholds in domainator hmmer_search were chosen for proteins, so scaling
    them by log(K)/log(20) keeps them at the same relative position on the score
    scale of a smaller alphabet.

    Returns 1.0 for the amino alphabet (leaving protein behavior untouched) and
    about 0.4628 for DNA and RNA.
    """
    return math.log(alphabet.K) / math.log(_AMINO_ALPHABET_K)


def _hmm_file_display_name(file) -> str:
    """Best-effort human readable name for something passed to pyhmmer.plan7.HMMFile."""
    if isinstance(file, (str, os.PathLike)):
        return os.fspath(file)
    return getattr(file, "name", repr(file))


class _DecompressingHMMFile(pyhmmer.plan7.HMMFile):
    """An HMMFile that also closes the decompressing handle it was built on.

    pyhmmer does not take ownership of a file object passed to HMMFile, and an
    HMMFile cannot carry an extra attribute, so the handle has to be held by a
    subclass to be closed deterministically.
    """

    def __init__(self, handle):
        self._handle = handle
        super().__init__(handle, db=False)

    def close(self):
        # Guard against closing twice: pyhmmer 0.12 aborts with a double free if
        # HMMFile.close() is called a second time. close() is also reached from
        # pyhmmer's __dealloc__, which can run during interpreter shutdown after
        # the attribute is already gone, hence the getattr.
        if not self.closed:
            super().close()
        handle = getattr(self, "_handle", None)
        if handle is not None:
            self._handle = None
            handle.close()


def open_hmm_file(file: Union[str, os.PathLike, IOBase]) -> pyhmmer.plan7.HMMFile:
    """Open an hmm file for reading, transparently handling gzip and BGZF.

    Two behaviors differ deliberately from handing the path straight to pyhmmer:

    Compression is detected from magic bytes rather than from the extension, and
    compressed files are read through gzip.open() rather than by path. HMMER's
    own dispatch is extension-based -- it pipes any path ending in ".gz" through
    "gzip -dc", so an uncompressed file merely named ".gz" fails outright, and it
    rejects BGZF entirely because of BGZF's FEXTRA header. gzip.open() reads both
    plain gzip and BGZF (BGZF is valid multi-member gzip) regardless of the name.

    Pressed hmmpress sidecars (.h3f/.h3i/.h3m/.h3p) are ignored: db=False. They
    are keyed to the .hmm only by filename, so an edited .hmm whose sidecars were
    not regenerated silently yields the *old* profiles with no warning. The load
    speedup they buy is negligible next to the profile-to-sequence comparisons
    that dominate a search, so the correct trade is to always read the .hmm.

    Args:
        file: path to, or open binary handle on, an hmm file. Handles are passed
            through as-is; compression and sidecar detection are both path-based.
            ``-`` or ``/dev/stdin`` reads standard input (uncompressed only).
    """
    if file in STDIN_PATHS:
        return pyhmmer.plan7.HMMFile(sys.stdin.buffer, db=False)
    if isinstance(file, (str, os.PathLike)):
        if detect_compression(file) is not None:
            return _DecompressingHMMFile(gzip.open(file, "rb"))
        if is_compressed_path(file):
            # Uncompressed content under a compressed name. Handing the path to
            # pyhmmer would make easel pipe it through "gzip -dc" and fail, so
            # open it ourselves and let the content decide.
            return _DecompressingHMMFile(open(file, "rb"))
    return pyhmmer.plan7.HMMFile(file, db=False)


def iter_hmm_names(file: Union[str, os.PathLike]) -> Iterable[str]:
    """Yield the profile names in an hmm file without parsing the profiles.

    Reads only the NAME lines, which is many times faster than building HMM
    objects and avoids holding the whole database in memory, for callers that
    want names only. A profile name cannot contain whitespace and NAME is a
    single line, so this yields exactly what hmm.name would.

    Handles gzip and BGZF, like open_hmm_file(). Raises no error on a file with
    no NAME lines; callers that need to distinguish "no profiles" from "not an
    hmm file after all" should check for an empty result.
    """
    opener = gzip.open if detect_compression(file) is not None else open
    with opener(file, "rb") as handle:
        for line in handle:
            if line.startswith(b"NAME"):
                parts = line.split(None, 1)
                if len(parts) == 2:
                    yield parts[1].strip().decode()


def iter_hmms(file: Union[str, os.PathLike, IOBase]) -> Iterable[pyhmmer.plan7.HMM]:
    """Iterate over the profiles in one hmm file, with a readable alphabet error.

    pyhmmer locks an HMMFile to the alphabet of its first profile and then raises a
    bare AlphabetMismatch, so a file holding more than one alphabet fails partway
    through with no indication of which file is at fault. This wrapper turns that
    into a message naming the file.
    """
    name = _hmm_file_display_name(file)
    try:
        with open_hmm_file(file) as hmm_file:
            yield from hmm_file
    except pyhmmer.errors.AlphabetMismatch as exc:
        raise _mixed_alphabet_in_file_error("input", name, exc) from None


def iter_hmms_with_alphabet(file: Union[str, os.PathLike, IOBase], role: str = "input") -> Tuple[Optional[pyhmmer.easel.Alphabet], Iterable[pyhmmer.plan7.HMM]]:
    """Open an hmm file and return (alphabet, iterator over all of its profiles).

    Reads the first profile eagerly so that its alphabet is known before the
    caller dispatches any work, then chains that profile back onto the rest of
    the file. This exists rather than calling peek_hmm_alphabet() first because
    the hmmer_* tools accept stdin, which cannot be reopened or rewound.

    Returns (None, empty iterator) for a file that opens but holds no profiles. Note
    that pyhmmer raises EOFError for a completely empty file rather than opening it.

    Args:
        file: path to, or open binary handle on, an hmm file.
        role: what the file is to the calling tool, interpolated into error messages.
    """
    name = _hmm_file_display_name(file)
    handle = open_hmm_file(file)
    try:
        first = handle.read()
    except pyhmmer.errors.AlphabetMismatch as exc:
        handle.close()
        raise _mixed_alphabet_in_file_error(role, name, exc) from None
    except BaseException:
        handle.close()
        raise
    if first is None:
        handle.close()
        return None, iter(())
    return first.alphabet, _iter_remaining_hmms(first, handle, role, name)


def _iter_remaining_hmms(first, handle, role, name):
    try:
        yield first
        while True:
            try:
                hmm = handle.read()
            except pyhmmer.errors.AlphabetMismatch as exc:
                raise _mixed_alphabet_in_file_error(role, name, exc, first.alphabet) from None
            if hmm is None:
                return
            yield hmm
    finally:
        handle.close()


def _mixed_alphabet_in_file_error(role, name, exc, first_alphabet=None) -> ValueError:
    seen = "" if first_alphabet is None else f" Its first profile is {alphabet_name(first_alphabet)}."
    return ValueError(
        f"The {role} hmm file '{name}' holds profiles with more than one alphabet.{seen} "
        f"HMMER locks an hmm file to the alphabet of its first profile, so such a file "
        f"cannot be read past that profile (pyhmmer reported: {exc}). "
        f"Split it by alphabet with domainator hmmer_select --alphabet."
    )


def peek_hmm_alphabet(file: Union[str, os.PathLike]) -> Optional[pyhmmer.easel.Alphabet]:
    """Return the alphabet of the first profile in an hmm file, or None if it has none.

    Only accepts paths: this reads the beginning of the file, so calling it on a
    stream would consume input that the caller still needs.
    """
    if not isinstance(file, (str, os.PathLike)):
        raise TypeError(
            "peek_hmm_alphabet only accepts file paths, because it consumes the start "
            f"of the file. Got: {type(file).__name__}."
        )
    with open_hmm_file(file) as hmm_file:
        for hmm in hmm_file:
            return hmm.alphabet
    return None


def common_hmm_alphabet(files: Iterable[Union[str, os.PathLike]], role: str = "input") -> Optional[pyhmmer.easel.Alphabet]:
    """Return the alphabet shared by a set of hmm files.

    Args:
        files: paths to hmm files. Streams are skipped, since they cannot be peeked
            at without consuming input the caller needs.
        role: what these files are to the calling tool ("input" or "reference"),
            interpolated into the error message so the user knows which argument to fix.

    Returns:
        The shared Alphabet, or None if none of the files held any profiles.

    Raises:
        ValueError: if two files use different alphabets. A single hmm file cannot
            hold more than one alphabet, so tools that concatenate profiles from
            several files would otherwise write output that cannot be read back.
    """
    alphabet = None
    alphabet_source = None
    for file in files:
        if not isinstance(file, (str, os.PathLike)):
            continue
        file_alphabet = peek_hmm_alphabet(file)
        if file_alphabet is None:
            continue
        if alphabet is None:
            alphabet = file_alphabet
            alphabet_source = os.fspath(file)
        elif file_alphabet != alphabet:
            raise ValueError(
                f"Mismatched alphabets among the {role} hmm files: "
                f"{alphabet_source} is {alphabet_name(alphabet)}, but "
                f"{os.fspath(file)} is {alphabet_name(file_alphabet)}. "
                "All profiles compared or written together must use the same alphabet."
            )
    return alphabet


def read_hmms(hmm_files:Iterable[Union[str,os.PathLike,IOBase]]) -> Dict[str, Dict[str,pyhmmer.plan7.HMM]]:
    """
        hmm_files: a list of paths to .hmm files

        returns:
            a dict of dicts of pyhmmer HMM objects
                db_name: hmm_name: HMM

    """
    out = dict()
    for file in hmm_files:
        if isinstance(file, IOBase):
            name = file.name
        else:
            name = file
        name = db_name_from_path(name)

        hmmer_models = OrderedDict() 
        for model in iter_hmms(file):
            model_name = model.name
            if model_name in hmmer_models:
                warnings.warn(f"multiple hmms with the same name ({model_name}) in file: {file}, only one will be used.")
            hmmer_models[model_name] = model
        if name in out:
            raise RuntimeError(f"Multiple hmm files with the same name, please combine the hmms into a single file, or rename one of the files. This is important to avoid searching the same domains twice, and for naming the source databases.")
        out[name] = hmmer_models
    return out


def read_pyhmmer_fastas(sequence_files):
    out = dict()
    for file_name in sequence_files:
        name = db_name_from_path(file_name)
        seqs_dict = dict()
        with pyhmmer.easel.SequenceFile(file_name, digital=True) as seq_file:
            for seq in seq_file:
                seq_name = seq.name
                if seq_name in seqs_dict:
                    warnings.warn(f"multiple reference sequences with the same name ({seq_name}) in file: {file_name}, only one will be used.")
                seqs_dict[seq_name] = seq
        if name in out:
            raise RuntimeError(f"Multiple reference sequence files with the same name, please combine the reference sequences into a single file, or rename one of the files. This is important to avoid searching the same domains twice, and for naming the source databases.")
        out[name] = seqs_dict
    return out


def read_pyhmmer_peptide_fastas(peptide_files):
    return read_pyhmmer_fastas(peptide_files)


def read_infernal_cms(cm_files):
    out = dict()
    for file_name in cm_files:
        name = db_name_from_path(file_name)
        cms_dict = OrderedDict()
        with pyinfernal.cm.CMFile(file_name) as cm_file:
            for cm in cm_file:
                cm_name = cm.name
                if cm_name in cms_dict:
                    warnings.warn(f"multiple covariance models with the same name ({cm_name}) in file: {file_name}, only one will be used.")
                cms_dict[cm_name] = cm
        if name in out:
            raise RuntimeError(f"Multiple cm files with the same name, please combine the models into a single file, or rename one of the files.")
        out[name] = cms_dict
    return out


def check_writable_hmm_path(path):
    """Raise if an hmm output path asks for a compression HMMER cannot read.

    Callers that stage output through a temporary file should call this on the
    real destination before doing any work, since the temporary name does not
    carry the destination's extension.
    """
    if path is None:
        return
    if Path(path).suffix.lower() in (".bgz", ".bgzf"):
        raise ValueError(
            f"Cannot write '{path}': BGZF-compressed hmm files are not supported, "
            "because HMMER cannot read them (easel rejects the BGZF header, so "
            "hmmsearch, hmmpress and Domainator itself would all fail to open the "
            "result). Use a '.gz' extension to write plain gzip instead."
        )


def open_writable_hmm_file(path, compressed=None):
    """Open a path for writing hmm profiles, using plain gzip for compressed outputs.

    Deliberately unlike open_writable_seqfile() above, which writes BGZF: HMMER
    cannot read a BGZF profile file. easel dispatches on the filename extension
    and pipes ".gz" paths through "gzip -dc", and it rejects BGZF's FEXTRA header
    outright, so a BGZF .hmm is unreadable by hmmsearch, by hmmpress, and by
    Domainator's own path-based reads. Plain gzip is read correctly by all of
    them. BGZF's advantage is seekability, and nothing seeks within a .hmm.

    ``compressed`` overrides the extension-based decision, for callers writing to
    a temporary path whose own extension does not reflect the destination.
    ``-`` and ``/dev/stdout`` return standard output (binary), which ``close()`` only flushes.
    """
    if path in STDOUT_PATHS:
        return _UnclosedStdout(sys.stdout.buffer)
    check_writable_hmm_path(path)
    if compressed is None:
        compressed = is_compressed_path(path)
    if compressed:
        return gzip.open(path, "wb")
    return open(path, "wb")
