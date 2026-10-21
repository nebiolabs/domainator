"""``domainator doctor``: which optional parts of Domainator this installation has.

Most of Domainator needs only its Python dependencies. A few commands also need an external
program or an optional package; this lists each one, the commands that need it, whether it
was found, and how to install it.
"""
import argparse
import importlib.util
import platform
import shutil
import sys
from typing import List, NamedTuple, Optional


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
    Requirement("tmbed", "find_features --tmbed", "pip install git+https://github.com/BernhoferM/TMbed.git"),
)


def _program_found(name: str) -> Optional[str]:
    return shutil.which(name)


def _package_found(name: str) -> Optional[str]:
    spec = importlib.util.find_spec(name)
    return None if spec is None else (spec.origin or "found")


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
    return lines


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(prog="domainator doctor", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.parse_args(argv)
    print("\n".join(report()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
