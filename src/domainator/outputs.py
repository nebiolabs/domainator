"""Writing a command's output files."""

import contextlib
import io
import os
from typing import Iterable, Optional, Sequence

from Bio import SeqIO
from domainator.output_guardrails import enforce_output_limit, make_temporary_output_path
from domainator.seqio import STDOUT_PATHS, is_compressed_path, open_writable_seqfile, write_genbank


class _Output:
    """One output path. A file is written to a temporary name and renamed into place on success."""

    def __init__(self, path: str):
        self.path = path
        if path in STDOUT_PATHS:
            self.temp_path = None
            self.handle = open_writable_seqfile(path)
        else:
            self.temp_path = make_temporary_output_path(path)
            self.handle = open_writable_seqfile(self.temp_path, compressed=is_compressed_path(path))

    def finish(self, success: bool):
        self.handle.close()
        if self.temp_path is None:
            return
        if success:
            os.replace(self.temp_path, self.path)
        elif os.path.exists(self.temp_path):
            os.remove(self.temp_path)


@contextlib.contextmanager
def output_file(path: str):
    """Open an output path for writing text, as a context manager.

    ``-`` and ``/dev/stdout`` write standard output. A file is written to a temporary name and
    renamed into place only if the block finishes without an exception, and a path ending in
    ``.gz`` or ``.bgz`` is compressed (BGZF, which gzip reads).
    """
    output = _Output(path)
    try:
        yield output.handle
    except BaseException:
        output.finish(success=False)
        raise
    output.finish(success=True)


class SeqOutputs:
    """Writes records to a GenBank output, a FASTA output, or both, in one pass.

    Use as a context manager: files appear at their paths only if the block finishes without
    an exception. With ``max_output_bytes``, raises ``OutputSizeLimitExceeded`` before writing a
    record that would take the outputs past the limit (``--max_output_gb``).

    Args:
        genbank: path of the GenBank output, ``-`` for standard output, or None.
        fasta: path of the FASTA output, ``-`` for standard output, or None.
        max_output_bytes: limit on the total bytes written to all outputs, or None.
        output_description, mitigation_options, extra_guidance: for the size-limit error message.
        default_molecule_type: GenBank molecule type for records that do not have one.
    """

    def __init__(self, genbank: Optional[str] = None, fasta: Optional[str] = None, *,
                 max_output_bytes: Optional[int] = None, output_description: Optional[str] = None,
                 mitigation_options: Sequence[str] = (), extra_guidance: Optional[str] = None,
                 default_molecule_type: str = "DNA"):
        self.genbank_path = genbank
        self.fasta_path = fasta
        self.max_output_bytes = max_output_bytes
        self.output_description = output_description or " and ".join(
            f"{kind} output '{path}'" for kind, path in (("GenBank", genbank), ("FASTA", fasta)) if path is not None)
        self.mitigation_options = list(mitigation_options)
        self.extra_guidance = extra_guidance
        self.default_molecule_type = default_molecule_type
        self.bytes_written = 0
        self._genbank = self._fasta = None

    def __enter__(self):
        if self.genbank_path is not None:
            self._genbank = _Output(self.genbank_path)
        if self.fasta_path is not None:
            self._fasta = _Output(self.fasta_path)
        return self

    def __exit__(self, exc_type, exc, traceback):
        for output in (self._genbank, self._fasta):
            if output is not None:
                output.finish(success=exc_type is None)
        return False

    def write(self, record) -> None:
        texts = []
        if self._genbank is not None:
            buffer = io.StringIO()
            write_genbank((record,), buffer, default_molecule_type=self.default_molecule_type)
            texts.append((self._genbank, buffer.getvalue()))
        if self._fasta is not None:
            buffer = io.StringIO()
            SeqIO.write(record, buffer, "fasta")
            texts.append((self._fasta, buffer.getvalue()))
        record_bytes = sum(len(text) for _, text in texts)  # GenBank and FASTA are ASCII
        enforce_output_limit(
            projected_bytes=self.bytes_written + record_bytes,
            max_output_bytes=self.max_output_bytes,
            output_description=self.output_description,
            mitigation_options=self.mitigation_options,
            extra_guidance=self.extra_guidance,
        )
        self.bytes_written += record_bytes
        for output, text in texts:
            output.handle.write(text)

    def write_all(self, records: Iterable) -> None:
        for record in records:
            self.write(record)


def write_seq_outputs(records: Iterable, genbank: Optional[str] = None, fasta: Optional[str] = None, **kwargs) -> int:
    """Write ``records`` to the given outputs (see :class:`SeqOutputs`) and return the bytes written."""
    with SeqOutputs(genbank, fasta, **kwargs) as outputs:
        outputs.write_all(records)
    return outputs.bytes_written
