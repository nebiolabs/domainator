[index](README.md)
# Domainator Basic Examples
These are simplest-case examples. For more complex examples, see the [companion repository](https://github.com/nebiolabs/domainator_examples).

## Download and prepare example data
These data are used in the examples below

```
wget https://zenodo.org/records/10989173/files/example_hmms.hmm # some various mostly unrelated hmm profiles
wget https://zenodo.org/records/10989173/files/ncbi_complete_subset.gb # about 800 Mb of various prokaryote genomes
wget https://zenodo.org/records/10989173/files/full_CDA_pfam_clan.hmm # HMM profiles for a specific family of related proteins
hmmer_select.py -i example_files/example_hmms.hmm --exact Sod_Fe_C -o Sod_Fe_C.hmm
```


## Run domainator on a single file

genbank input with CDSs already annotated
```
domainate.py --cpu 4 -i example_files/pDONR201_genemark.gb -r example_files/example_hmms.hmm -o annotated.gb --max_overlap 0.6
plot_contigs.py -i annotated.gb --html contigs_plot.html
```

protein fasta input
```
domainate.py --cpu 4 -i example_files/example_peptides.fasta -r example_files/example_hmms.hmm -o annotated.gb --max_overlap 0.6
plot_contigs.py -i annotated.gb --html contigs_plot.html
```

nucleotide fasta input
```
domainate.py --gene_call all --cpu 4 -i example_files/pDONR201_multi.fasta -r example_files/example_hmms.hmm -o annotated.gb --max_overlap 0.6
plot_contigs.py -i annotated.gb --html contigs_plot.html
```

fasta reference (runs phmmer instead of hmmscan)
```
domainate.py --cpu 4 -i example_files/pDONR201_genemark.gb -r example_files/example_peptides.fasta -o annotated.gb --max_overlap 0.6
plot_contigs.py -i annotated.gb --html contigs_plot.html
```



## Run domainator on multiple input files

You can have multiple input files. 
In this example, one of the inputs is a nucleotide genbank file with CDS annotations already present and the other is a fasta file with no annotations. setting `--gene_call unannotated` will find CDSs in any contig that doesn't have at least one CDS annotation already.
```
domainate.py --gene_call unannotated --cpu 4 -o annotated.gb -r example_files/example_hmms.hmm --max_overlap 0.6 -i example_files/pDONR201_genemark.gb example_files/pDONR201_multi.fasta
```

You can also search multiple domain databases at once by supplying multiple hmm or protein fasta files via the `-r` option. Not all of the reference sequences have to be of the same type. Some can be protein sequences and some can be hmm files in the same call to domainate.

You can also run domainate sequentially, meaning you can run an already annotated genbank file through the domainator program again with a different domain database. Options such as `--max_overlap` and `--max_domains` are only applied to newly added annotations.

## Visualizing output
```
summary_report.py -i annotated.gb --html domain_summary.html
enum_report.py -i annotated.gb --by contig --name --taxname superkingdom genus species self --length --domains --html enum_report.html -o enum_report.tsv
color_genbank.py -i annotated.gb -o annotated_colored.gb --color_domains
plot_contigs.py -i annotated_colored.gb --html contigs_plot.html
```

## Searching a database of genomes and returning genome neighborhoods

References (queries) can be either protein hmm profiles or protein fasta files.
```
domain_search.py --max_hits 100 --cpu 4 -i ncbi_complete_subset.gb -r example_files/ -e 1e-10 -o neighborhoods.gb --cds_range 10
plot_contigs.py -i neighborhoods.gb --html neighborhoods_plot.html
```



## Making a sequences similarity network
```
seq_dist.py -i proteins.fasta -r proteins.fasta --dense_text efi_scores.tsv --mode efi_score
enum_report.py -i proteins.fasta --length --sequence -o metadata.tsv
build_ssn.py -i efi_scores.tsv --metadata metadata.tsv --lb 50 --xgmml cytoscape_file.lb50.xgmml --cluster_tsv clusters.tsv 
```
In practice, you should probably use `--sparse efi_scores.hdf` for output, instead of `--dense_text`. Particularly if you have a lot of sequences. Otherwise the distance matrix file will be huge. Text output is nice for readability and compatibility with other programs though.

For the SSN, you will probably need to try several `--lb` settings to find one that gives a convenient number of clusters.

## Parallelizing domainate or domain_search over multiple nodes of a compute cluster

`partition_seqfile.py` followed by `domainate.py` or `domain_search.py` using the `--offset` and `--recs_to_read` options.