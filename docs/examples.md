[index](README.md)
# Domainator Basic Examples
These are simplest-case examples. For more complex examples, see the [companion repository](https://github.com/nebiolabs/domainator_examples).

Many of the examples below are also included in an executable [Google Colab notebook](https://colab.research.google.com/github/nebiolabs/domainator_examples/blob/main/colab_notebooks/Domainator.ipynb)

## Download and prepare example data
These data are used in the examples below

```bash
wget https://zenodo.org/records/10989173/files/ncbi_complete_subset.gb # about 800 Mb of various prokaryote genomes
wget https://zenodo.org/records/10999789/files/example_files.tar.gz
tar -xvf example_files.tar.gz
hmmer_select.py -i example_files/example_hmms.hmm --exact Sod_Fe_C -o Sod_Fe_C.hmm
```


## Run domainator on a single file

genbank input with CDSs already annotated
```bash
domainate.py --cpu 4 -i example_files/pDONR201_genemark.gb -r example_files/example_hmms.hmm -o annotated.gb --max_overlap 0.6
plot_contigs.py -i annotated.gb --html contigs_plot.html
```

protein fasta input
```bash
domainate.py --cpu 4 -i example_files/example_peptides.fasta -r example_files/example_hmms.hmm -o annotated.gb --max_overlap 0.6
plot_contigs.py -i annotated.gb --html contigs_plot.html
```

nucleotide fasta input
```bash
domainate.py  --fasta_type nucleotide --gene_call all --cpu 4 -i example_files/pDONR201_multi.fasta -r example_files/example_hmms.hmm -o annotated.gb --max_overlap 0.6
plot_contigs.py -i annotated.gb --html contigs_plot.html
```

fasta reference (runs phmmer instead of hmmscan)
```bash
domainate.py --cpu 4 -i example_files/pDONR201_genemark.gb -r example_files/example_peptides.fasta -o annotated.gb --max_overlap 0.6
plot_contigs.py -i annotated.gb --html contigs_plot.html
```



## Run domainator on multiple input files

You can have multiple input files. 
In this example, one of the inputs is a nucleotide genbank file with CDS annotations already present and the other is a fasta file with no annotations. setting `--gene_call unannotated` will find CDSs in any contig that doesn't have at least one CDS annotation already.
```bash
domainate.py --fasta_type nucleotide --gene_call unannotated --cpu 4 -o annotated.gb -r example_files/example_hmms.hmm --max_overlap 0.6 -i example_files/pDONR201_genemark.gb example_files/pDONR201_multi.fasta
```

You can also search multiple domain databases at once by supplying multiple hmm or protein fasta files via the `-r` option. Not all of the reference sequences have to be of the same type. Some can be protein sequences and some can be hmm files in the same call to domainate.

You can also run domainate sequentially, meaning you can run an already annotated genbank file through the domainator program again with a different domain database. Options such as `--max_overlap` and `--max_domains` are only applied to newly added annotations.

## Visualizing output
```bash
summary_report.py -i annotated.gb --html domain_summary.html
enum_report.py -i annotated.gb --by contig --definition --taxname superkingdom genus species self --length --domains --html enum_report.html -o enum_report.tsv
color_genbank.py -i annotated.gb -o annotated_colored.gb --color_domains
plot_contigs.py -i annotated_colored.gb --html contigs_plot.html
```

## Searching a database of genomes

Input can be CDS-annotated genbank files, protein genbank files, protein fasta files, or even nucleotide fasta files (if `--fasta_type nucleotide` and `--gene_call` are specified). 

References (queries) can be either protein hmm profiles or protein fasta files.
```bash
# find hits and extract their neighborhoods
domain_search.py --max_hits 100 --cpu 4 -i ncbi_complete_subset.gb -r Sod_Fe_C.hmm -e 1e-10 -o neighborhoods.gb --cds_range 10

# annotate the contigs with additional domain annotations other than the query
domainate.py -i neighborhoods.gb -o annotated_neighborhoods.gb -r example_files/example_hmms.hmm --max_overlap 0.6

# generate various kinds of reports of the results
summary_report.py -i annotated_neighborhoods.gb -o /dev/stdout --html summary_report.html --taxonomy
enum_report.py -i annotated_neighborhoods.gb -o /dev/stdout --html enum_report.html --definition --length --taxname superkingdom genus self --domains --domain_descriptions --domain_search
color_genbank.py -i annotated_neighborhoods.gb -o annotated_neighborhoods_colored.gb --color_domains
plot_contigs.py -i annotated_neighborhoods_colored.gb --html neighborhoods_plot.html 
```

Output protein sequences
```bash
domain_search.py --max_hits 100 --translate --cpu 4 -i ncbi_complete_subset.gb -r Sod_Fe_C.hmm -e 1e-10 -o protein_hits.gb 
domainate.py -i protein_hits.gb  -o annotated_proteins.gb -r example_files/example_hmms.hmm --max_overlap 0.6

# generate various kinds of reports
summary_report.py -i annotated_proteins.gb -o /dev/stdout --html summary_report.html --taxonomy
enum_report.py -i annotated_proteins.gb -o /dev/stdout --html enum_report.html --definition --length --taxname superkingdom genus self --domains --domain_descriptions --domain_search
color_genbank.py -i annotated_proteins.gb -o annotated_proteins_colored.gb --color_domains
plot_contigs.py -i annotated_proteins_colored.gb --html proteins_plot.html 
```

## Making a sequence similarity network or tree
```bash
enum_report.py -i example_files/FeSOD_20.gb --length --sequence -o metadata.tsv

# --dense_text gives a nice human readable score matrix files, but it can be very large, frequently you'd prefer --sparse or --dense for more efficient files.
# --mode score can be later converted to various other kinds of scores, like efi_score and efi_score_dist, etc. If you only want one output, you can set --mode to that directly.
seq_dist.py -i example_files/FeSOD_20.gb -r example_files/FeSOD_20.gb --sparse scores.sparse.hdf5 --mode score
transform_matrix.py -i scores.sparse.hdf5 --dense efi_score_dist.dense.hdf5 --mode efi_score_dist

# Domainator doesn't have a newick renderer, but the xgmml output can be opened in cytoscape, you'll need to apply a layout after opening it.
build_tree.py -i efi_score_dist.dense.hdf5 --newick tree.nwk --xgmml tree.xgmml --metadata metadata.tsv


transform_matrix.py -i scores.sparse.hdf5 --sparse efi_scores.sparse.hdf5 --mode efi_score

# shows a score histogram to give you guess at where to set the --lb for the ssn. You'll probably still have to iterate a few times on build_ssn.py
# --include_mst_knn adds projected MST_KNN edge counts, but it uses extra memory for neighbor rankings and threshold-by-k counts. Use --profile_stages to inspect those sizes for large matrices.
matrix_report.py -i efi_scores.sparse.hdf5 -o /dev/stdout

# you can open sequence_similarity_network.xgmml in Cytoscape, and apply force directed layout to get a nice visualization
# --cluster creates a new metadata column called 'SSN_cluster', following connectivity. 
# --color_by can be any column from the metadata table, or SSN_cluster (if --cluster is set)
# --color_table_out saves the color table to a tsv file.
# To thin out a hairball, sparsify: --mst keeps just the maximum spanning tree, --knn K keeps
# each node's K nearest neighbors, and --mst --knn K (or the shorthand --mst_knn K) keeps both.
# Only the options that include --mst preserve the clusters of the full graph; --knn on its own
# can break one cluster into several.
build_ssn.py -i efi_scores.sparse.hdf5 --xgmml sequence_similarity_network.xgmml --lb 50 --cluster --cluster_tsv clusters.tsv --color_by SSN_cluster --color_table_out clusters_colors.tsv --metadata metadata.tsv

# generates a nice figure legend from the cluster labels.
color_table_to_legend.py -i clusters_colors.tsv --svg cluster_legend.svg
```
## Making a neighborhood similarity network or tree

Same as above, except use `compare_contigs.py` instead of `seq_dist.py`, this will compare contigs based on Jaccard index and/or adjacency index of their annotated domain contents.

## Building and comparing nucleotide profiles

The `hmmer_*` programs work on DNA and RNA profiles as well as amino acid ones. Because
HMMER locks an `.hmm` file to the alphabet of its first profile, every profile in one file
must share an alphabet.

```bash
# Build a DNA profile from a nucleotide alignment. The alphabet is inferred when it can be,
# but pass --alphabet to be sure: an alignment of only A/C/G/T is also valid protein.
# --window_length sets the profile's MAXL field, which nhmmer uses to window long targets.
hmmer_build.py -i dna_alignment.afa --alphabet dna --name my_dna_profile --window_length 500 -o my_dna_profile.hmm

# Confirm what you got. --alphabet reports amino, DNA, or RNA.
hmmer_report.py -i my_dna_profile.hmm --alphabet --length --consensus -o dna_profile_metadata.tsv

# Compare nucleotide profiles to each other. Scores are on the nucleotide scale, which runs
# at roughly half the amino acid scale, so protein-tuned thresholds are about twice as
# strict here. hmmer_search.py scales its default --score_cutoff for you.
hmmer_compare.py -i my_dna_profiles.hmm -r my_dna_profiles.hmm -o dna_scores.tsv --alignments
hmmer_search.py -i my_dna_profiles.hmm -r reference_dna_profiles.hmm -o dna_hits.hmm

# Split a mixed set of profile files by alphabet.
hmmer_select.py -i *.hmm --alphabet dna -o dna_only.hmm
hmmer_select.py -i *.hmm --alphabet amino -o amino_only.hmm
```

A nucleotide `.hmm` built this way can be used directly as a `domainate.py` or
`domain_search.py` reference, where it is searched against whole nucleotide contigs with
nhmmer rather than against individual CDSs.

## Making a profile tree or profile similarity network
```bash
# write some metadata to a tab-separated file
hmmer_report.py -i example_files/full_CDA_pfam_clan.hmm --length --acc --desc --consensus -o CDS_metadata.tsv

# --mode score can be later converted to various other kinds of scores, like efi_score and efi_score_dist, etc. If you only want one output, you can set --mode to that directly.
seq_dist.py -i example_files/full_CDA_pfam_clan.hmm -r example_files/full_CDA_pfam_clan.hmm --algorithm hmmer_compare --mode score --sparse CDA_score.sparse.hdf5
transform_matrix.py -i CDA_score.sparse.hdf5 --dense CDA_efi_score_dist.dense.hdf5 --mode efi_score_dist

# Domainator doesn't have a newick renderer, but the xgmml output can be opened in cytoscape, you'll need to apply a layout after opening it.
build_tree.py -i CDA_efi_score_dist.dense.hdf5 --newick tree.nwk --xgmml tree.xgmml --metadata CDS_metadata.tsv


transform_matrix.py -i CDA_score.sparse.hdf5 --sparse CDA_efi_score.sparse.hdf5 --mode efi_score
# shows a score histogram to give you guess at where to set the --lb for the ssn. You'll probably still have to iterate a few times on build_ssn.py
# --include_mst_knn adds projected MST_KNN edge counts, but it uses extra memory for neighbor rankings and threshold-by-k counts. Use --profile_stages to inspect those sizes for large matrices.
matrix_report.py -i CDA_efi_score.sparse.hdf5 -o /dev/stdout

# you can open profile_similarity_network.xgmml in Cytoscape, and apply force directed layout to get a nice visualization
# --cluster creates a new metadata column called 'SSN_cluster', following connectivity. 
# --color_by can be any column from the metadata table, or SSN_cluster (if --cluster is set)
# --color_table_out saves the color table to a tsv file.
build_ssn.py -i CDA_efi_score.sparse.hdf5 --xgmml profile_similarity_network.xgmml --lb 10 --cluster --cluster_tsv clusters.tsv --color_by SSN_cluster --color_table_out clusters_colors.tsv --metadata CDS_metadata.tsv

# generates a nice figure legend from the cluster labels.
color_table_to_legend.py -i clusters_colors.tsv --svg cluster_legend.svg
```

## Slicing and subsetting genbank files

 - `select_by_cds.py` is useful for extracting gene neighborhoods or individual genes from larger annotated contigs, using advanced selection criteria.
 - `select_by_contig.py` is useful for extracting whole contigs from larger files, using advanced selection criteria.
 - `extract_domains.py`, `extract_peptides.py`, `extract_unannotated.py`, and `deduplicate_genbank.py` are also useful for their specialized functions.

 Check their documentation, using `-h`. And also see the advanced workflows at the [companion repository](https://github.com/nebiolabs/domainator_examples).


## Parallelizing domainate or domain_search over multiple nodes of a compute cluster

`partition_seqfile.py` followed by `domainate.py` or `domain_search.py` using the `--offset` and `--recs_to_read` options.
