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
domainator hmmer_select -i example_files/example_hmms.hmm --exact Sod_Fe_C -o Sod_Fe_C.hmm
```


## Run domainator on a single file

genbank input with CDSs already annotated
```bash
domainator domainate --cpu 4 -i example_files/pDONR201_genemark.gb -r example_files/example_hmms.hmm -o annotated.gb --max_overlap 0.6
domainator plot_contigs -i annotated.gb --html_out contigs_plot.html
```

protein fasta input
```bash
domainator domainate --cpu 4 -i example_files/example_peptides.fasta -r example_files/example_hmms.hmm -o annotated.gb --max_overlap 0.6
domainator plot_contigs -i annotated.gb --html_out contigs_plot.html
```

nucleotide fasta input
```bash
domainator domainate  --alphabet nucleotide --gene_call all --cpu 4 -i example_files/pDONR201_multi.fasta -r example_files/example_hmms.hmm -o annotated.gb --max_overlap 0.6
domainator plot_contigs -i annotated.gb --html_out contigs_plot.html
```

fasta reference (runs phmmer instead of hmmscan)
```bash
domainator domainate --cpu 4 -i example_files/pDONR201_genemark.gb -r example_files/example_peptides.fasta -o annotated.gb --max_overlap 0.6
domainator plot_contigs -i annotated.gb --html_out contigs_plot.html
```



## Run domainator on multiple input files

You can have multiple input files. 
In this example, one of the inputs is a nucleotide genbank file with CDS annotations already present and the other is a fasta file with no annotations. setting `--gene_call unannotated` will find CDSs in any contig that doesn't have at least one CDS annotation already.
```bash
domainator domainate --alphabet nucleotide --gene_call unannotated --cpu 4 -o annotated.gb -r example_files/example_hmms.hmm --max_overlap 0.6 -i example_files/pDONR201_genemark.gb example_files/pDONR201_multi.fasta
```

You can also search multiple domain databases at once by supplying multiple hmm or protein fasta files via the `-r` option. Not all of the reference sequences have to be of the same type. Some can be protein sequences and some can be hmm files in the same call to domainate.

You can also run domainate sequentially, meaning you can run an already annotated genbank file through the domainator program again with a different domain database. Options such as `--max_overlap` and `--max_hits_per_cds` are only applied to newly added annotations.

## Visualizing output
```bash
domainator summary_report -i annotated.gb --html_out domain_summary.html
domainator enum_report -i annotated.gb --by contig --columns definition taxname:superkingdom:genus:species:self length domains --html_out enum_report.html -o enum_report.tsv
domainator color_genbank -i annotated.gb -o annotated_colored.gb --target domains
domainator plot_contigs -i annotated_colored.gb --html_out contigs_plot.html
```

## Searching a database of genomes

Input can be CDS-annotated genbank files, protein genbank files, protein fasta files, or even nucleotide fasta files (if `--alphabet nucleotide` and `--gene_call` are specified). 

References (queries) can be either protein hmm profiles or protein fasta files.
```bash
# find hits and extract their neighborhoods
domainator domain_search --max_hits 100 --cpu 4 -i ncbi_complete_subset.gb -r Sod_Fe_C.hmm -e 1e-10 -o neighborhoods.gb --cds_range 10

# annotate the contigs with additional domain annotations other than the query
domainator domainate -i neighborhoods.gb -o annotated_neighborhoods.gb -r example_files/example_hmms.hmm --max_overlap 0.6

# generate various kinds of reports of the results
domainator summary_report -i annotated_neighborhoods.gb -o /dev/stdout --html_out summary_report.html --taxonomy
domainator enum_report -i annotated_neighborhoods.gb -o /dev/stdout --html_out enum_report.html --columns definition length taxname:superkingdom:genus:self domains domain_descriptions domain_search
domainator color_genbank -i annotated_neighborhoods.gb -o annotated_neighborhoods_colored.gb --target domains
domainator plot_contigs -i annotated_neighborhoods_colored.gb --html_out neighborhoods_plot.html 
```

Output protein sequences
```bash
domainator domain_search --max_hits 100 --translate --cpu 4 -i ncbi_complete_subset.gb -r Sod_Fe_C.hmm -e 1e-10 -o protein_hits.gb 
domainator domainate -i protein_hits.gb  -o annotated_proteins.gb -r example_files/example_hmms.hmm --max_overlap 0.6

# generate various kinds of reports
domainator summary_report -i annotated_proteins.gb -o /dev/stdout --html_out summary_report.html --taxonomy
domainator enum_report -i annotated_proteins.gb -o /dev/stdout --html_out enum_report.html --columns definition length taxname:superkingdom:genus:self domains domain_descriptions domain_search
domainator color_genbank -i annotated_proteins.gb -o annotated_proteins_colored.gb --target domains
domainator plot_contigs -i annotated_proteins_colored.gb --html_out proteins_plot.html 
```

## Making a sequence similarity network or tree
```bash
domainator enum_report -i example_files/FeSOD_20.gb --columns length sequence -o metadata.tsv

# --dense_text_out gives a nice human readable score matrix files, but it can be very large, frequently you'd prefer --hdf_out for more efficient files.
# --mode score can be later converted to various other kinds of scores, like efi_score and efi_score_dist, etc. If you only want one output, you can set --mode to that directly.
domainator seq_dist -i example_files/FeSOD_20.gb -r example_files/FeSOD_20.gb --hdf_out scores.hdf5 --mode score
domainator transform_matrix -i scores.hdf5 --hdf_out efi_score_dist.hdf5 --mode efi_score_dist

# Domainator doesn't have a newick renderer, but the xgmml output can be opened in cytoscape, you'll need to apply a layout after opening it.
domainator build_tree -i efi_score_dist.hdf5 --newick_out tree.nwk --xgmml_out tree.xgmml --metadata metadata.tsv


domainator transform_matrix -i scores.hdf5 --hdf_out efi_scores.hdf5 --mode efi_score

# shows a score histogram to give you guess at where to set the --threshold for the ssn. You'll probably still have to iterate a few times on domainator build_ssn
# --knn_projection adds the edge counts domainator build_ssn --mst --knn K would keep at each threshold, but it uses extra memory for neighbor rankings and threshold-by-k counts. Use --profile_stages to inspect those sizes for large matrices.
domainator matrix_report -i efi_scores.hdf5 -o /dev/stdout

# you can open sequence_similarity_network.xgmml in Cytoscape, and apply force directed layout to get a nice visualization
# --cluster creates a new metadata column called 'SSN_cluster', following connectivity. 
# --color_by can be any column from the metadata table, or SSN_cluster (if --cluster is set)
# --color_table_out saves the color table to a tsv file.
# To thin out a hairball, sparsify: --mst keeps just the maximum spanning tree, --knn K keeps
# each node's K nearest neighbors, and --mst --knn K keeps both. domainator seq_dist, domainator compare_contigs,
# domainator structure_dist, domainator kmer_dist, and domainator transform_matrix take the same two options.
# Only the options that include --mst preserve the clusters of the full graph; --knn on its own
# can break one cluster into several.
domainator build_ssn -i efi_scores.sparse.hdf5 --xgmml_out sequence_similarity_network.xgmml --threshold 50 --cluster --cluster_tsv_out clusters.tsv --color_by SSN_cluster --color_table_out clusters_colors.tsv --metadata metadata.tsv

# generates a nice figure legend from the cluster labels.
domainator color_table_to_legend -i clusters_colors.tsv --svg_out cluster_legend.svg
```
## Making a neighborhood similarity network or tree

Same as above, except use `domainator compare_contigs` instead of `domainator seq_dist`, this will compare contigs based on Jaccard index and/or adjacency index of their annotated domain contents.

## Building and comparing nucleotide profiles

The `hmmer_*` programs work on DNA and RNA profiles as well as amino acid ones. Because
HMMER locks an `.hmm` file to the alphabet of its first profile, every profile in one file
must share an alphabet.

```bash
# Build a DNA profile from a nucleotide alignment. The alphabet is inferred when it can be,
# but pass --alphabet to be sure: an alignment of only A/C/G/T is also valid protein.
# --window_length sets the profile's MAXL field, which nhmmer uses to window long targets.
domainator hmmer_build -i dna_alignment.afa --alphabet dna --name my_dna_profile --window_length 500 -o my_dna_profile.hmm

# Confirm what you got. --alphabet reports amino, DNA, or RNA.
domainator hmmer_report -i my_dna_profile.hmm --columns alphabet length consensus -o dna_profile_metadata.tsv

# Compare nucleotide profiles to each other. Scores are on the nucleotide scale, which runs
# at roughly half the amino acid scale, so protein-tuned thresholds are about twice as
# strict here. domainator hmmer_search scales its default --min_score for you.
domainator hmmer_compare -i my_dna_profiles.hmm -r my_dna_profiles.hmm -o dna_scores.tsv --alignments_out dna_alignments.txt
domainator hmmer_search -i my_dna_profiles.hmm -r reference_dna_profiles.hmm -o dna_hits.hmm

# Split a mixed set of profile files by alphabet.
domainator hmmer_select -i *.hmm --alphabet dna -o dna_only.hmm
domainator hmmer_select -i *.hmm --alphabet amino -o amino_only.hmm
```

A nucleotide `.hmm` built this way can be used directly as a `domainator domainate` or
`domainator domain_search` reference, where it is searched against whole nucleotide contigs with
nhmmer rather than against individual CDSs.

## Making a profile tree or profile similarity network
```bash
# write some metadata to a tab-separated file
domainator hmmer_report -i example_files/full_CDA_pfam_clan.hmm --columns length acc desc consensus -o CDS_metadata.tsv

# --mode score can be later converted to various other kinds of scores, like efi_score and efi_score_dist, etc. If you only want one output, you can set --mode to that directly.
domainator seq_dist -i example_files/full_CDA_pfam_clan.hmm -r example_files/full_CDA_pfam_clan.hmm --algorithm hmmer_compare --mode score --hdf_out CDA_score.hdf5
domainator transform_matrix -i CDA_score.hdf5 --hdf_out CDA_efi_score_dist.hdf5 --mode efi_score_dist

# Domainator doesn't have a newick renderer, but the xgmml output can be opened in cytoscape, you'll need to apply a layout after opening it.
domainator build_tree -i CDA_efi_score_dist.hdf5 --newick_out tree.nwk --xgmml_out tree.xgmml --metadata CDS_metadata.tsv


domainator transform_matrix -i CDA_score.hdf5 --hdf_out CDA_efi_score.hdf5 --mode efi_score
# shows a score histogram to give you guess at where to set the --threshold for the ssn. You'll probably still have to iterate a few times on domainator build_ssn
# --knn_projection adds the edge counts domainator build_ssn --mst --knn K would keep at each threshold, but it uses extra memory for neighbor rankings and threshold-by-k counts. Use --profile_stages to inspect those sizes for large matrices.
domainator matrix_report -i CDA_efi_score.hdf5 -o /dev/stdout

# you can open profile_similarity_network.xgmml in Cytoscape, and apply force directed layout to get a nice visualization
# --cluster creates a new metadata column called 'SSN_cluster', following connectivity. 
# --color_by can be any column from the metadata table, or SSN_cluster (if --cluster is set)
# --color_table_out saves the color table to a tsv file.
domainator build_ssn -i CDA_efi_score.sparse.hdf5 --xgmml_out profile_similarity_network.xgmml --threshold 10 --cluster --cluster_tsv_out clusters.tsv --color_by SSN_cluster --color_table_out clusters_colors.tsv --metadata CDS_metadata.tsv

# generates a nice figure legend from the cluster labels.
domainator color_table_to_legend -i clusters_colors.tsv --svg_out cluster_legend.svg
```

## Slicing and subsetting genbank files

 - `domainator select_by_cds` is useful for extracting gene neighborhoods or individual genes from larger annotated contigs, using advanced selection criteria.
 - `domainator select_by_contig` is useful for extracting whole contigs from larger files, using advanced selection criteria.
 - `domainator rotate_contigs` sets the origin of circular contigs, such as plasmids or phages, relative to a domain, so that related contigs start at the same place and in the same orientation. For example, `domainator rotate_contigs -i plasmids.gb -r CcdB.hmm --anchor start -o rotated.gb`.
 - `domainator extract_domains`, `domainator extract_peptides`, `domainator extract_unannotated`, and `domainator deduplicate_genbank` are also useful for their specialized functions.

 Check their documentation, using `-h`. And also see the advanced workflows at the [companion repository](https://github.com/nebiolabs/domainator_examples).


## Parallelizing domainate or domain_search over multiple nodes of a compute cluster

Split the database into shards and index them with `domainator format_db --shards N --index`, then run one `domainator domainate` or `domainator domain_search` job per shard. Within a single machine, `domainator domain_search` already parallelizes over partitions of each input file (`--cpu`), and an index lets it skip the offset scan.

## More commands, one example each

These assume `annotated.gb` is the output of `domainator domainate` on nucleotide contigs,
and `annotated_proteins.gb` the same for proteins. Every command's options are listed in the
[command-line reference](cli/README.md).

Select, filter and edit records:

```bash
# Keep only the contigs that have both a CcdA and a CcdB domain, and no CAT domain
domainator select_by_contig -i annotated.gb --domain_expr "CcdA & CcdB & ~CAT" -o selected.gb

# Keep domain hits with E-values of 1e-10 or less that overlap a better hit by at most 50%
domainator filter_domains -i annotated.gb -e 1e-10 --max_overlap 0.5 -o filtered.gb

# Trim 5 kb from both ends of each contig
domainator trim_contigs -i annotated.gb --kb_range 5 -o trimmed.gb

# Trim CDSs that have a CcdB domain off the ends of each contig
domainator trim_contigs -i annotated.gb --domains CcdB -o trimmed.gb

# Extract the regions listed in a table of contig, start, end (1-based and inclusive), and optionally strand
domainator select_by_coord -i annotated.gb --coords regions.tsv -o regions.gb

# Write the 10 best domain_search hits, best first (for domain_search output written without --max_hits)
domainator sort_contigs -i hits.gb --max_hits 10 -o top_hits.gb

# Make sequence names safe for Newick trees and drop sequences with non-canonical residues, logging what changed
domainator clean_sequences -i proteins.fasta --alphabet protein --clean_name --filter_by_aa -o clean.gb --log_out clean.log
```

Extract sequences:

```bash
# The proteins with a CcdB domain, as a protein GenBank file and a FASTA file
domainator extract_peptides -i annotated.gb --domains CcdB -o ccdb_proteins.gb --fasta_out ccdb_proteins.fasta

# The parts of each protein at least 50 residues long that no domain covers (linkers, unknown domains)
domainator extract_unannotated -i annotated_proteins.gb --min_length 50 -o unannotated.gb

# Every record as FASTA
domainator genbank_to_fasta -i annotated.gb --fasta_out annotated.fasta
```

Annotate features other than domains:

```bash
# Mark every match to a PROSITE pattern (here two cysteines two residues apart)
domainator find_features -i annotated_proteins.gb --motif "CX(2)C" -o motifs.gb
```

Compare and project:

```bash
# k-mer similarity between every pair of proteins
domainator kmer_dist -i proteins.fasta -r proteins.fasta -o kmer_scores.hdf5

# Similarity of contigs by their domain content, and the contigs sorted by it
domainator compare_contigs -i annotated.gb --ji 0.5 --ai 0.5 -o contig_scores.hdf5 --genbank_out sorted_contigs.gb

# A 2D projection of a similarity matrix, as a Cytoscape network and an interactive plot (which saves PNG and SVG pictures)
domainator build_projection -i scores.hdf5 --algorithm pca -o projection.xgmml --html_out projection.html
```

Build a covariance model:

```bash
# From a Stockholm alignment with a consensus secondary structure (experimental; needs Infernal)
domainator infernal_build -i my_rna.sto --name my_rna -o my_rna.cm
```
