[index](README.md)
# Domainator Basic Examples
These are simplest-case examples. For more complex examples, see the [companion repository](https://github.com/nebiolabs/domainator_examples).

## Run domainator on a single file
```
domainator.py --cpu 4 -i your_genbank.gb -r Pfam-A.hmm -o domainator_output.gb --max_overlap 0.6
```

## Run domainator on multiple input files
```
domainator.py --cpu 4 -o domainator_output.gb -r Pfam-A.hmm --max_overlap 0.6 -i folder_with_your_genbanks/*.gb
```
or
```
domainator.py --cpu 4 -o domainator_output.gb -r Pfam-A.hmm --max_overlap 0.6 -i your_genbank_1.gb your_genbank_2.gb 
```

You can also search multiple domain databases at once by supplying multiple hmm files via the `-r` option

You can also run the domainator program sequentially, meaning you can run an already annotated genbank file through the domainator program again with a different domain database.

## Visualizing output
```
summary_report.py -i domainator_output.gb --html domain_summary.html
enum_report.py -i domainator_output.gb --by contig --name --taxname superkingdom genus species self --length --domains --html enum_report.html -o enum_report.tsv
plot_contigs.py -i domainator_output.gb --html contigs_plot.html
```

## Searching a database of genomes and returning genome neighborhoods
```
domain_search.py --max_hits 100 --cpu 4 -i genome_database.gb -r query_profile.hmm -e 1e-10 -o neighborhoods.gb --cds_range 10
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