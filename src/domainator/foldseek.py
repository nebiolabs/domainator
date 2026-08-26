try:
    from esmologs.ESM2_to_3Di import ESM2_to_3Di
    from esmologs.predict_from_ESM2_to_3Di import convert_batch
    from esmologs.fasta2foldseek import fasta2foldseek
    import torch
except ImportError: 
    pass 

import os
import shutil
import tempfile
import subprocess
from typing import List, Iterable, Tuple, Union, Iterator
from collections import namedtuple

# define a named tuple for hits with fields "query,target,qheader,theader,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
FoldseekHit = namedtuple("Hit", ["query","target","qheader","theader","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits", "qlen", "tlen"])

MAX_PROTEIN_SIZE = 2500


def resolve_foldseek_path(foldseek_path=None) -> str:
    """Resolve the foldseek executable, raising a descriptive error if it is missing."""
    if foldseek_path is not None:
        return foldseek_path
    resolved = shutil.which("foldseek")
    if resolved is None:
        raise RuntimeError(
            "Could not find 'foldseek' on PATH. Install foldseek "
            "(for example 'conda install -c conda-forge -c bioconda foldseek') "
            "or pass an explicit path."
        )
    return resolved


def _run(argv, description, env=None):
    """Run a foldseek subcommand, raising RuntimeError with stderr on failure.

    Uses subprocess.run rather than Popen+wait: the previous form read stderr only
    after the process exited, which deadlocks once stderr fills the pipe buffer.
    """
    completed = subprocess.run(argv, capture_output=True, env=env)
    if completed.returncode != 0:
        raise RuntimeError(
            f"{description} exited with code {completed.returncode}:\n"
            + completed.stderr.decode("utf-8", errors="replace")
        )
    return completed

def search(database_path, proteins, foldseek, cpu, E, device=None, foldseek_path=None) -> Iterable[FoldseekHit]:
    """
    Args:
        foldseek_path: explicit path to the foldseek executable. If None, it is
            looked up on PATH.
        device: which device foldseek should run the search on. Either None or "cpu" for
            CPU-only search, or a CUDA device string ("cuda", "cuda:0", "cuda:1", ...) to
            enable GPU-accelerated search. For a GPU search the target database must have
            been built in a GPU-compatible (padded) format (see `foldseek makepaddedseqdb`).
    """
    foldseek_bin = resolve_foldseek_path(foldseek_path)
    with tempfile.TemporaryDirectory() as tmpdirname:
        out_base_name = tmpdirname + "/output"
        protein_fasta_name = tmpdirname + "/protein.fasta"
        threedi_fasta_name = tmpdirname + "/threedi.fasta"
        foldseek_tmpfolder = tmpdirname + "/foldseek_tmpfolder"
        aln_path = tmpdirname + "/aln"

        num_seqs = 0
        with open(protein_fasta_name, "w") as protein_f:
            with open(threedi_fasta_name, "w") as threedi_f:
                for i, foldseek_seq in enumerate(foldseek):
                    if foldseek_seq is None:
                        continue
                    num_seqs += 1

                    protein = proteins[i]
                    protein = protein.textize()

                    protein_f.write(f">{protein.name} {protein.description}\n{protein.sequence}\n")
                    threedi_f.write(foldseek_seq + "\n")
        
        if num_seqs == 0:
            return # no sequences to search, yield nothing
        
        fasta2foldseek(protein_fasta_name, threedi_fasta_name, out_base_name)
        foldseek_options = [foldseek_bin, "search", out_base_name, database_path, aln_path, foldseek_tmpfolder, "-e", str(E)]
        if cpu > 0 and cpu is not None:
            foldseek_options += ["--threads", str(cpu)]
        foldseek_env = None
        if device is not None and device != "cpu":
            # GPU-accelerated search. foldseek selects the GPU via the CUDA_VISIBLE_DEVICES
            # environment variable rather than a command line flag, so map a "cuda:N" device
            # string onto that variable.
            foldseek_options += ["--gpu", "1"]
            if ":" in device:
                foldseek_env = os.environ.copy()
                foldseek_env["CUDA_VISIBLE_DEVICES"] = device.split(":", 1)[1]
        _run(foldseek_options, "foldseek", env=foldseek_env)

        _run(
            [foldseek_bin, "convertalis", "--format-output",
             "query,target,qheader,theader,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen",
             out_base_name, database_path, aln_path, foldseek_tmpfolder + "/results.tsv"],
            "convertalis",
        )

        with open(foldseek_tmpfolder + "/results.tsv", "r") as f:
            for line in f:
                yield FoldseekHit(*line.strip().split("\t"))


class foldseekBuilder():
    def __init__(self, device="cuda:0", checkpoint=None):
        self.device = device
        self.checkpoint = checkpoint
        self.model = ESM2_to_3Di("esm2_t36_3B_UR50D", torch.load(checkpoint, map_location=device))
        self.checkpoint=checkpoint
        self.device=device
        self.model.to(self.device)
        self.model.eval()

    def __call__(self, name:str, prot:str) -> bytes:
        if len(prot) > MAX_PROTEIN_SIZE: #TODO: maybe warn?
            return None
        # skip if contains non-amino acid characters
        if prot.strip("ACDEFGHIKLMNPQRSTVWY") != "": #TODO: maybe warn?
            return None
        predicted_seqs = convert_batch(self.model, [prot], device=self.device)
        return f">{name}\n{predicted_seqs[0]}"