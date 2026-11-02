"""``domainator doctor``: which optional parts of Domainator this installation has.

Most of Domainator needs only its Python dependencies. A few commands also need an external
program or an optional package; this lists each one, the commands that need it, whether it
was found, and how to install it. It also reports whether PyTorch can use an NVIDIA GPU.
"""
import argparse
import importlib.util
import platform
import re
import shutil
import subprocess
import sys
from typing import List, NamedTuple, Optional, Tuple


class Requirement(NamedTuple):
    name: str
    needed_by: str
    install: str


EXTERNAL_PROGRAMS = (
    Requirement("diamond", "seq_dist --algorithm diamond_* (the default), deduplicate_genbank --algorithm diamond",
                "conda install -c bioconda diamond"),
    Requirement("cd-hit", "deduplicate_genbank --algorithm cd-hit (the default; cd-hit-est for nucleotides)",
                "conda install -c bioconda cd-hit"),
    Requirement("usearch", "deduplicate_genbank --algorithm usearch", "conda install -c bioconda usearch"),
    Requirement("foldseek", "domainate --foldseek, structure_* --aligner foldseek (the default)",
                "conda install -c conda-forge -c bioconda foldseek"),
    Requirement("reseek", "structure_* --aligner reseek", "build it from https://github.com/rcedgar/reseek"),
    Requirement("cmbuild", "infernal_build", "conda install -c bioconda infernal"),
    Requirement("datasets", "db_download of NCBI genomes (--download_backend datasets, the default)",
                "conda install -c conda-forge ncbi-datasets-cli"),
    Requirement("unzip", "db_download of NCBI genomes", "your system's package manager"),
)

OPTIONAL_PACKAGES = (
    Requirement("umap", "build_projection --algorithm umap", 'pip install "domainator[projection]"'),
    Requirement("torch", "domainate --foldseek on sequences (ESM-2), find_features --tmbed",
                'pip install "domainator[structure]"'),
    Requirement("esmologs", "domainate --foldseek on sequences (ESM-2)", "see docs/esm_3b_foldseek.md"),
    Requirement("tmbed", "find_features --tmbed", 'see "Optional: TMbed" in README.md (needs Python 3.14 or older)'),
)

GPU_NEEDED_BY = ("domainate --esm2_3di_device, --prostt5_device and --foldseek_device (a GPU build of foldseek); "
                 "find_features --device (TMbed)")
PYTORCH_INSTALL = "https://pytorch.org/get-started/locally/"


def _program_found(name: str) -> Optional[str]:
    return shutil.which(name)


def _package_found(name: str) -> Optional[str]:
    spec = importlib.util.find_spec(name)
    return None if spec is None else (spec.origin or "found")


def _nvidia_gpus() -> Optional[Tuple[List[str], Optional[str]]]:
    """The GPUs nvidia-smi lists, and the newest CUDA version their driver supports. None without nvidia-smi."""
    nvidia_smi = shutil.which("nvidia-smi")
    if nvidia_smi is None:
        return None
    try:
        listing = subprocess.run([nvidia_smi, "-L"], capture_output=True, text=True, timeout=30).stdout
        summary = subprocess.run([nvidia_smi], capture_output=True, text=True, timeout=30).stdout
    except (OSError, subprocess.SubprocessError):
        return [], None
    gpus = [line.split(" (UUID")[0] for line in listing.splitlines() if line.startswith("GPU ")]
    cuda = re.search(r"CUDA Version:\s*([\d.]+)", summary)
    return gpus, cuda.group(1) if cuda else None


def _torch_cuda() -> Optional[Tuple[str, Optional[str], int]]:
    """torch's version, the CUDA version it was built for (None for a CPU-only build), and how many
    GPUs it can use. None when torch is not installed."""
    if _package_found("torch") is None:
        return None
    import torch
    return torch.__version__, torch.version.cuda, torch.cuda.device_count() if torch.cuda.is_available() else 0


def _cuda_version(version: str) -> Tuple[int, ...]:
    return tuple(int(part) for part in version.split("."))


def _torch_status(driver_cuda: Optional[str]) -> str:
    try:
        torch_cuda = _torch_cuda()
    except Exception as error:  # a broken torch install
        return f"cannot be imported ({error})"
    if torch_cuda is None:
        return "not installed"
    version, built_cuda, usable = torch_cuda
    if built_cuda is None:
        return f"{version} is a CPU-only build; for a GPU, install a CUDA build ({PYTORCH_INSTALL})"
    if usable:
        return f"{version}, built for CUDA {built_cuda}, can use {usable} GPU{'' if usable == 1 else 's'}"
    if driver_cuda and _cuda_version(built_cuda) > _cuda_version(driver_cuda):
        return (f"{version} is built for CUDA {built_cuda}, newer than the driver supports, so it runs on the CPU. "
                f"Install a build for CUDA {driver_cuda} or older ({PYTORCH_INSTALL}), or update the NVIDIA driver.")
    return f"{version}, built for CUDA {built_cuda}, cannot use a GPU"


def gpu_report() -> List[str]:
    gpus = _nvidia_gpus()
    driver_cuda = None
    if gpus is None:
        found = "none found (nvidia-smi is not on PATH)"
    elif not gpus[0]:
        found = "none found (nvidia-smi lists no GPUs)"
    else:
        names, driver_cuda = gpus
        found = ", ".join(names) + (f"; the driver supports CUDA {driver_cuda} and older" if driver_cuda else "")
    return ["", "GPU:",
            f"  NVIDIA GPU: {found}",
            f"  torch:      {_torch_status(driver_cuda)}",
            f"  needed by:  {GPU_NEEDED_BY}"]


def report() -> List[str]:
    from domainator import __version__
    from domainator.seqio import native_parser_available

    lines = [f"Domainator {__version__}, Python {platform.python_version()} on {platform.platform()}", ""]
    lines.append("Native (Rust) parser: " + ("installed" if native_parser_available() else
                 "NOT installed; GenBank and FASTA parsing falls back to Biopython (slower). "
                 "See the README's 'Native (Rust) acceleration' section."))
    for heading, requirements, found in (("External programs", EXTERNAL_PROGRAMS, _program_found),
                                          ("Optional packages", OPTIONAL_PACKAGES, _package_found)):
        lines += ["", f"{heading}:"]
        for requirement in requirements:
            where = found(requirement.name)
            status = f"found ({where})" if where else f"missing; install: {requirement.install}"
            lines.append(f"  {requirement.name:<9} {status}")
            lines.append(f"  {'':<9} needed by: {requirement.needed_by}")
    return lines + gpu_report()


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(prog="domainator doctor", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.parse_args(argv)
    print("\n".join(report()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
