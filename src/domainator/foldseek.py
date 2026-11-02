import os
import tempfile
from typing import Iterable
from collections import namedtuple
from domainator.external import import_optional, resolve_binary, run_command

# define a named tuple for hits with fields "query,target,qheader,theader,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
FoldseekHit = namedtuple("Hit", ["query","target","qheader","theader","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits", "qlen", "tlen"])

MAX_PROTEIN_SIZE = 2500
STANDARD_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def predictable(prot: str) -> bool:
    """Whether a 3Di predictor is run on ``prot``: at most MAX_PROTEIN_SIZE residues, all standard."""
    return len(prot) <= MAX_PROTEIN_SIZE and prot.strip(STANDARD_AMINO_ACIDS) == ""


def _gpu_options(device):
    """foldseek's options and environment for running on ``device``: "cpu" (or None), or a CUDA
    device string ("cuda", "cuda:0", ...). foldseek selects the GPU via the CUDA_VISIBLE_DEVICES
    environment variable rather than a command line flag, so "cuda:N" is mapped onto that."""
    if device is None or device == "cpu":
        return [], None
    env = None
    if ":" in device:
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = device.split(":", 1)[1]
    return ["--gpu", "1"], env


def _import_optional(module: str):
    """Import one of the structure-prediction packages (esmologs, torch), which are optional."""
    return import_optional(module, "See docs/esm_3b_foldseek.md for how to install the structure-search dependencies.")


def search(database_path, proteins, foldseek, cpu, E, device=None, foldseek_path=None,
           prostt5_weights=None, prostt5_device=None) -> Iterable[FoldseekHit]:
    """
    Args:
        foldseek: per protein, its predicted 3Di sequence (as a fasta record), or None to skip
            it. Under ProstT5 an entry only marks a protein for the search (see ProstT5Builder).
        prostt5_weights: foldseek's ProstT5 weights (from 'foldseek databases ProstT5'). If
            given, foldseek predicts the 3Di sequences itself, with ProstT5.
        prostt5_device: the device ProstT5 runs on: "cpu", or a CUDA device string, as for
            ``device``.
        foldseek_path: explicit path to the foldseek executable. If None, it is
            looked up on PATH.
        device: which device foldseek should run the search on. Either None or "cpu" for
            CPU-only search, or a CUDA device string ("cuda", "cuda:0", "cuda:1", ...) to
            enable GPU-accelerated search. For a GPU search the target database must have
            been built in a GPU-compatible (padded) format (see `foldseek makepaddedseqdb`).
    """
    foldseek_bin = resolve_binary("foldseek", foldseek_path, install_hint="Install foldseek (for example 'conda install -c conda-forge -c bioconda foldseek').")
    with tempfile.TemporaryDirectory() as tmpdirname:
        out_base_name = tmpdirname + "/output"
        protein_fasta_name = tmpdirname + "/protein.fasta"
        threedi_fasta_name = tmpdirname + "/threedi.fasta"
        foldseek_tmpfolder = tmpdirname + "/foldseek_tmpfolder"
        aln_path = tmpdirname + "/aln"

        threads = ["--threads", str(cpu)] if cpu is not None and cpu > 0 else []
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
                    if prostt5_weights is None:
                        threedi_f.write(foldseek_seq + "\n")
        
        if num_seqs == 0:
            return # no sequences to search, yield nothing
        
        if prostt5_weights is None:
            _import_optional("esmologs.fasta2foldseek").fasta2foldseek(protein_fasta_name, threedi_fasta_name, out_base_name)
        else:
            gpu_options, gpu_env = _gpu_options(prostt5_device)
            run_command([foldseek_bin, "createdb", protein_fasta_name, out_base_name,
                         "--prostt5-model", str(prostt5_weights)] + threads + gpu_options,
                        "foldseek createdb", env=gpu_env)
        gpu_options, foldseek_env = _gpu_options(device)
        foldseek_options = [foldseek_bin, "search", out_base_name, database_path, aln_path, foldseek_tmpfolder,
                            "-e", str(E)] + threads + gpu_options
        run_command(foldseek_options, "foldseek search", env=foldseek_env)

        run_command(
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
        torch = _import_optional("torch")
        ESM2_to_3Di = _import_optional("esmologs.ESM2_to_3Di").ESM2_to_3Di
        self._convert_batch = _import_optional("esmologs.predict_from_ESM2_to_3Di").convert_batch
        self.device = device
        self.checkpoint = checkpoint
        self.model = ESM2_to_3Di("esm2_t36_3B_UR50D", torch.load(checkpoint, map_location=device))
        self.model.to(self.device)
        self.model.eval()

    def __call__(self, name:str, prot:str) -> bytes:
        if not predictable(prot): #TODO: maybe warn?
            return None
        predicted_seqs = self._convert_batch(self.model, [prot], device=self.device)
        return f">{name}\n{predicted_seqs[0]}"

class ProstT5Builder():
    """The ProstT5 counterpart of foldseekBuilder. foldseek itself predicts the 3Di sequences,
    in ``search``, so this only selects the proteins to search, by the same rule."""

    def __call__(self, name: str, prot: str):
        return "" if predictable(prot) else None
