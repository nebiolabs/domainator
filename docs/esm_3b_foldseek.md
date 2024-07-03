Domainator can support sequence searches at high sensitivity across deep evolutionary distances by leveraging the ESM-2 3B 3Di model, described in the manuscript:
[https://www.biorxiv.org/content/10.1101/2023.07.26.550718v1](https://www.biorxiv.org/content/10.1101/2023.07.26.550718v1)

# Installation
In addition to Domainator, the esmologs package ([https://github.com/seanrjohnson/esmologs](https://github.com/seanrjohnson/esmologs)) and pytorch with CUDA support must be installed. The easiest way to accomplish this is to create a new conda environment with esmologs, and install domainator into that environment.

You also need to download the ESM-2 3B 3Di fine tuning checkpoint, from [https://zenodo.org/record/8174960](https://zenodo.org/record/8174960)

## download the checkpoint

```bash
wget https://zenodo.org/record/8174960/files/ESM-2_3B_3Di.pt
```

## Install via conda

```bash
git clone https://github.com/seanrjohnson/esmologs.git
cd esmologs

conda env create --name domainator_esmologs --file conda_env.yml

cd ..
git clone https://github.com/nebiolabs/domainator.git
cd domainator
conda env update --name domainator_esmologs --file conda_env.yml

conda activate domainator_esmologs
pytest test
cd ..
```

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

# Using ESM-2 3B 3Di with domainate.py

In this workflow, we first create a reference Foldseek 3Di database, and then use domainator to annotate contigs from that database

(Note that domain_search.py with ESM-2 3B 3Di is not yet supported)

## conda
```bash
conda activate domainator_esmologs

# convert a reference file to 3di
predict_from_ESM2_to_3Di.py -i domainator/test/data/foldseek/FeSOD_20.fasta -o FeSOD_20.3di.fasta --weights ESM-2_3B_3Di.pt --device cuda:0

# convert the amino acid and 3di fasta files into a foldseek database
fasta2foldseek.py --aa domainator/test/data/foldseek/FeSOD_20.fasta --tdi FeSOD_20.3di.fasta -o FeSOD

# run domainate.py with the foldseek reference database. In this example, our query is the same file we used to make the database, but it could be any fasta or genbank file.
domainate.py -i domainator/test/data/foldseek/FeSOD_20.fasta -o FeSOD_all_to_all_3Di.gb --foldseek FeSOD --esm2_3Di_weights ESM-2_3B_3Di.pt --esm2_3Di_device cuda:0
```

## Apptainer/Singularity

```bash

# convert a reference file to 3di
apptainer exec --nv domainator/domainator_esmologs.sif predict_from_ESM2_to_3Di.py -i domainator/test/data/foldseek/FeSOD_20.fasta -o FeSOD_20.3di.fasta --weights ESM-2_3B_3Di.pt --device cuda:0

# convert the amino acid and 3di fasta files into a foldseek database
apptainer exec domainator/domainator_esmologs.sif fasta2foldseek.py --aa domainator/test/data/foldseek/FeSOD_20.fasta --tdi FeSOD_20.3di.fasta -o FeSOD

# run domainate.py with the foldseek reference database. In this example, our query is the same file we used to make the database, but it could be any fasta or genbank file.
apptainer exec --nv domainator/domainator_esmologs.sif domainate.py -i domainator/test/data/foldseek/FeSOD_20.fasta -o FeSOD_all_to_all_3Di.gb --foldseek FeSOD --esm2_3Di_weights ESM-2_3B_3Di.pt --esm2_3Di_device cuda:0
```

