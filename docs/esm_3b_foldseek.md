[index](README.md)
# ESM-2 3B 3Di and foldseek integration

Domainator can support sequence searches at high sensitivity across deep evolutionary distances by leveraging the ESM-2 3B 3Di model, described in the manuscript:
[https://www.biorxiv.org/content/10.1101/2023.07.26.550718v1](https://www.biorxiv.org/content/10.1101/2023.07.26.550718v1)

This page covers annotating *sequences* that have no structure, by predicting their 3Di
alphabet with ESM-2 3B. If you already have structures (pdb/cif files, or a foldseek
database built from them), use the structure tools described in
[Structure search and annotation](structure_search.md) instead: they align real coordinates,
support TM-align and the TM-score/lDDT/RMSD metrics, and need only the foldseek binary — no
PyTorch, no esmologs, and no checkpoint.

# ProstT5: the same workflow without PyTorch

foldseek (release 10 or later) can predict 3Di itself, with the ProstT5 model, running inside
the foldseek binary. Download its weights once, and pass them to `domainate` in place of the
ESM-2 3B 3Di checkpoint:

```bash
foldseek databases ProstT5 prostt5_weights tmp

domainator domainate -i domainator/test/data/foldseek/FeSOD_20.fasta -o FeSOD_all_to_all_ProstT5.gb --foldseek domainator/test/data/foldseek/FeSOD --prostt5_weights prostt5_weights
```

`--prostt5_device cuda:0` runs ProstT5 on a GPU, which needs a GPU build of foldseek. Proteins
are selected by the same rule as with ESM-2 3B 3Di (at most 2,500 residues, standard amino
acids only). ESM-2 3B 3Di and ProstT5 3Di give similar search sensitivity (see the
[eLife paper](https://elifesciences.org/articles/91415), Fig. 4), but they are different
models, so for the best scores build the reference database with the same predictor as the
queries (`foldseek createdb refs.fasta refs --prostt5-model prostt5_weights`), or from
structures. On `FeSOD_20` against the ESM-2 3B 3Di `FeSOD` database, ProstT5 queries find
the same 400 hits as ESM-2 3B 3Di queries, with bit scores about 10% lower.

# Installation
In addition to Domainator, the esmologs package ([https://github.com/seanrjohnson/esmologs](https://github.com/seanrjohnson/esmologs)), fair-esm and PyTorch with CUDA support must be installed. They can be installed into the `domainator` conda environment from the README's Installation section.

You also need to download the ESM-2 3B 3Di fine tuning checkpoint, from [https://zenodo.org/record/8174960](https://zenodo.org/record/8174960)

## download the checkpoint

```bash
wget https://zenodo.org/record/8174960/files/ESM-2_3B_3Di.pt
```

## Install into the domainator conda environment

```bash
conda activate domainator

# PyTorch: choose the build for your NVIDIA driver. The "CUDA Version" at the top of nvidia-smi's
# output is the newest CUDA the driver supports; install a build for that version or older (the
# builds are listed at https://pytorch.org/get-started/locally/). For example, for CUDA 12.4:
pip install torch --index-url https://download.pytorch.org/whl/cu124

pip install "fair-esm>=2.0" "esmologs @ git+https://github.com/seanrjohnson/esmologs.git@0e937758ebb68983a3049b64e62307e524c10ac4"

domainator doctor   # torch and esmologs should be found, and under GPU, torch should be able to use the GPU
```

pip's default PyTorch is built for a recent CUDA version (CUDA 13 for recent releases). With an
older driver it can't use the GPU, and `--esm2_3di_device cuda:0` (the default) fails;
`domainator doctor` says when this is the problem.

## install via Apptainer/Singularity

```bash
git clone https://github.com/nebiolabs/domainator.git
cd domainator

apptainer build domainator_esmologs.sif domainator_esmologs.def

# if using wsl, you need to use --nvccli. In other linux, --nv also works. These flags make the GPU visible to the container.
apptainer shell --nvccli domainator_esmologs.sif

# in the apptainer shell
cd /opt/domainator
pytest test
exit # or ctrl + d 
```

# Using ESM-2 3B 3Di with domainator domainate

In this workflow, we first create a reference Foldseek 3Di database, and then use domainator to annotate contigs from that database

(Note that domainator domain_search with ESM-2 3B 3Di is not supported. To search a large structure database with a small number of queries, use [domainator structure_search](structure_search.md).)

## conda
```bash
conda activate domainator

# convert a reference file to 3di
predict_from_ESM2_to_3Di.py -i domainator/test/data/foldseek/FeSOD_20.fasta -o FeSOD_20.3di.fasta --weights ESM-2_3B_3Di.pt --device cuda:0

# convert the amino acid and 3di fasta files into a foldseek database
fasta2foldseek.py --aa domainator/test/data/foldseek/FeSOD_20.fasta --tdi FeSOD_20.3di.fasta -o FeSOD

# run domainator domainate with the foldseek reference database. In this example, our query is the same file we used to make the database, but it could be any fasta or genbank file.
domainator domainate -i domainator/test/data/foldseek/FeSOD_20.fasta -o FeSOD_all_to_all_3Di.gb --foldseek FeSOD --esm2_3di_weights ESM-2_3B_3Di.pt --esm2_3di_device cuda:0
```

## Apptainer/Singularity

```bash

# convert a reference file to 3di
apptainer exec --nv domainator/domainator_esmologs.sif predict_from_ESM2_to_3Di.py -i domainator/test/data/foldseek/FeSOD_20.fasta -o FeSOD_20.3di.fasta --weights ESM-2_3B_3Di.pt --device cuda:0

# convert the amino acid and 3di fasta files into a foldseek database
apptainer exec domainator/domainator_esmologs.sif fasta2foldseek.py --aa domainator/test/data/foldseek/FeSOD_20.fasta --tdi FeSOD_20.3di.fasta -o FeSOD

# run domainator domainate with the foldseek reference database. In this example, our query is the same file we used to make the database, but it could be any fasta or genbank file.
apptainer exec --nv domainator/domainator_esmologs.sif domainator domainate -i domainator/test/data/foldseek/FeSOD_20.fasta -o FeSOD_all_to_all_3Di.gb --foldseek FeSOD --esm2_3di_weights ESM-2_3B_3Di.pt --esm2_3di_device cuda:0
```

