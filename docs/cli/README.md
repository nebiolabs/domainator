[index](../README.md)
# Command-line reference

Generated from the commands' parsers by `scripts/docs/cli_reference.py`. `domainator <command> --help` shows the same options.

## Annotate and search sequences and structures

- [domainate](domainate.md): Annotate contigs or proteins with hits to HMM, sequence, or CM references.
- [domain_search](domain_search.md): Search a large database with a few references and extract the hits.
- [find_features](find_features.md): Annotate transmembrane regions, signal peptides, and other sequence features.
- [structure_domainate](structure_domainate.md): Annotate protein structures with hits to reference structures.
- [structure_search](structure_search.md): Search a structure database with a few reference structures.

## Select, filter, and edit annotated records

- [clean_sequences](clean_sequences.md): Clean sequence names and residues, and filter out sequences with non-canonical residues.
- [color_genbank](color_genbank.md): Add color qualifiers to CDS and domain features.
- [deduplicate_genbank](deduplicate_genbank.md): Reduce redundancy by clustering sequences (cd-hit, usearch, diamond, or exact hashes).
- [extract_domains](extract_domains.md): Extract annotated domains as their own records.
- [extract_peptides](extract_peptides.md): Extract translated CDSs as protein records, keeping their domain annotations.
- [extract_unannotated](extract_unannotated.md): Extract the regions of proteins not covered by any domain annotation.
- [filter_domains](filter_domains.md): Filter domain annotations by e-value, overlap, database, or name.
- [rotate_contigs](rotate_contigs.md): Set the origin of circular contigs relative to a domain.
- [select_by_cds](select_by_cds.md): Extract contigs or neighborhoods around selected CDSs.
- [select_by_contig](select_by_contig.md): Select whole contigs by their domain content and other properties.
- [select_by_coord](select_by_coord.md): Extract contig regions listed in a coordinate table.
- [sort_contigs](sort_contigs.md): Sort contigs by their Domain_Search score.
- [trim_contigs](trim_contigs.md): Trim the ends of contigs.

## Convert between formats

- [genbank_to_fasta](genbank_to_fasta.md): Convert GenBank files to FASTA.
- [structure_to_genbank](structure_to_genbank.md): Convert protein structures into protein GenBank records.

## Build, select, and compare HMM and CM profiles

- [hmmer_build](hmmer_build.md): Build an HMM profile from a multiple sequence alignment.
- [infernal_build](infernal_build.md): Build a covariance model (CM) from a multiple sequence alignment (experimental).
- [hmmer_select](hmmer_select.md): Extract a subset of profiles from HMM files.
- [hmmer_search](hmmer_search.md): Select HMM profiles by their similarity to reference profiles.
- [hmmer_compare](hmmer_compare.md): Score and align pairs of HMM profiles.

## Compute similarity and distance matrices

- [seq_dist](seq_dist.md): Compute similarity scores between sequences or profiles.
- [kmer_dist](kmer_dist.md): Compute exact k-mer similarity scores between sequence records.
- [structure_dist](structure_dist.md): Compute similarity or distance matrices between protein structures.
- [compare_contigs](compare_contigs.md): Compute similarities between gene neighborhoods from their domain content and order.
- [transform_matrix](transform_matrix.md): Convert score matrices into normalized scores or distances.

## Build, explore, and split similarity networks, trees, and projections

- [build_ssn](build_ssn.md): Build a sequence similarity network as a Cytoscape file.
- [build_ssn_viewer](build_ssn_viewer.md): Build a Domainator Similarity Network Viewer bundle and HTML page.
- [ssn_navigator](ssn_navigator.md): Query a similarity network bundle: clusters, members, and metadata.
- [split_by_cluster](split_by_cluster.md): Split a sequence or profile file into one file per cluster of a similarity matrix.
- [build_tree](build_tree.md): Build a UPGMA tree from a distance or similarity matrix.
- [build_projection](build_projection.md): Project a matrix into 2D or 3D with UMAP, t-SNE, or PCA.

## Tables and summaries

- [enum_report](enum_report.md): Tabulate properties of each contig, CDS, or domain.
- [hmmer_report](hmmer_report.md): Tabulate properties of each HMM profile.
- [summary_report](summary_report.md): Summarize the domain and taxonomy content of a dataset.
- [matrix_report](matrix_report.md): Summarize a similarity matrix.

## Figures

- [plot_contigs](plot_contigs.md): Draw contigs and their annotations as an interactive HTML page.
- [color_table_to_legend](color_table_to_legend.md): Draw an SVG legend from a color table.

## Download and prepare databases

- [db_download](db_download.md): Download sequence databases in formats Domainator reads.
- [format_db](format_db.md): Shard, compress, and index a database for domain_search.
