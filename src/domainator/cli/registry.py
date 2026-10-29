"""The static catalog of Domainator commands.

This is the single source of the command list for the ``domainator`` dispatcher, its
help text, and shell completion. It is plain data so that listing or dispatching
commands imports no tool module; ``domainator <command>`` imports exactly one.
"""

from typing import Dict, NamedTuple, Tuple


class Tool(NamedTuple):
    module: str
    """Module under the ``domainator`` package that defines ``build_parser()`` and ``main(argv)``."""
    category: str
    summary: str


CATEGORIES: Tuple[Tuple[str, str], ...] = (
    ("annotate", "Annotate and search sequences and structures"),
    ("edit", "Select, filter, and edit annotated records"),
    ("convert", "Convert between formats"),
    ("profile", "Build, select, and compare HMM and CM profiles"),
    ("compare", "Compute similarity and distance matrices"),
    ("network", "Build, explore, and split similarity networks, trees, and projections"),
    ("report", "Tables and summaries"),
    ("plot", "Figures"),
    ("database", "Download and prepare databases"),
)
"""(key, heading) pairs, in the order ``domainator --help`` lists them."""


TOOLS: Dict[str, Tool] = {
    # annotate
    "domainate": Tool("domainate", "annotate", "Annotate contigs or proteins with reference hits, retaining records without hits."),
    "domain_search": Tool("domain_search", "annotate", "Select matching sequence records or hit neighborhoods; optionally translate CDS hits."),
    "find_features": Tool("find_features", "annotate", "Annotate transmembrane regions, signal peptides, and other sequence features."),
    "structure_domainate": Tool("structure_domainate", "annotate", "Annotate structure-derived protein records, retaining records without hits."),
    "structure_search": Tool("structure_search", "annotate", "Select matching structure-derived protein records and mark their best-hit regions."),
    # edit
    "clean_sequences": Tool("clean_sequences", "edit", "Clean sequence names and residues, and filter out sequences with non-canonical residues."),
    "color_genbank": Tool("color_genbank", "edit", "Add color qualifiers to CDS and domain features."),
    "deduplicate_genbank": Tool("deduplicate_genbank", "edit", "Reduce redundancy by clustering sequences (cd-hit, usearch, diamond, or exact hashes)."),
    "extract_domains": Tool("extract_domains", "edit", "Extract annotated domains as their own records."),
    "extract_peptides": Tool("extract_peptides", "edit", "Extract translated CDSs as protein records, keeping their domain annotations."),
    "extract_unannotated": Tool("extract_unannotated", "edit", "Extract the regions of proteins not covered by any domain annotation."),
    "filter_domains": Tool("filter_domains", "edit", "Filter domain annotations by e-value, overlap, database, or name."),
    "rotate_contigs": Tool("rotate_contigs", "edit", "Set the origin of circular contigs relative to a domain."),
    "select_by_cds": Tool("select_by_cds", "edit", "Extract contigs or neighborhoods around selected CDSs."),
    "select_by_contig": Tool("select_by_contig", "edit", "Select whole contigs by their domain content and other properties."),
    "select_by_coord": Tool("select_by_coord", "edit", "Extract contig regions listed in a coordinate table."),
    "sort_contigs": Tool("sort_contigs", "edit", "Sort contigs by their Domain_Search score."),
    "trim_contigs": Tool("trim_contigs", "edit", "Trim the ends of contigs."),
    # convert
    "genbank_to_fasta": Tool("genbank_to_fasta", "convert", "Convert GenBank files to FASTA."),
    "structure_to_genbank": Tool("structure_to_genbank", "convert", "Convert protein structures into protein GenBank records."),
    # profile
    "hmmer_build": Tool("hmmer_build", "profile", "Build an HMM profile from a multiple sequence alignment."),
    "infernal_build": Tool("infernal_build", "profile", "Build a covariance model (CM) from a multiple sequence alignment (experimental)."),
    "hmmer_select": Tool("hmmer_select", "profile", "Extract a subset of profiles from HMM files."),
    "hmmer_search": Tool("hmmer_search", "profile", "Select HMM profiles by their similarity to reference profiles."),
    "hmmer_compare": Tool("hmmer_compare", "profile", "Score and align pairs of HMM profiles."),
    # compare
    "seq_dist": Tool("seq_dist", "compare", "Compute similarity scores between sequences or profiles."),
    "kmer_dist": Tool("kmer_dist", "compare", "Compute exact k-mer similarity scores between sequence records."),
    "structure_dist": Tool("structure_dist", "compare", "Compute similarity or distance matrices between protein structures."),
    "compare_contigs": Tool("compare_contigs", "compare", "Compute similarities between gene neighborhoods from their domain content and order."),
    "transform_matrix": Tool("transform_matrix", "compare", "Convert score matrices into normalized scores or distances."),
    # network
    "build_ssn": Tool("build_ssn", "network", "Build a sequence similarity network for Cytoscape or the Domainator Viewer."),
    "ssn_navigator": Tool("ssn_navigator", "network", "Query a similarity network bundle: clusters, members, and metadata."),
    "split_by_cluster": Tool("split_by_cluster", "network", "Split a sequence or profile file into one file per cluster of a similarity matrix."),
    "build_tree": Tool("build_tree", "network", "Build a UPGMA tree from a distance or similarity matrix."),
    "build_projection": Tool("build_projection", "network", "Project a matrix into 2D or 3D with UMAP, t-SNE, or PCA."),
    "viewer": Tool("viewer", "network", "Write the Domainator Viewer, empty or opening a network, tree, or projection bundle."),
    # report
    "enum_report": Tool("enum_report", "report", "Tabulate properties of each contig, CDS, or domain."),
    "hmmer_report": Tool("hmmer_report", "report", "Tabulate properties of each HMM profile."),
    "summary_report": Tool("summary_report", "report", "Summarize the domain and taxonomy content of a dataset."),
    "matrix_report": Tool("matrix_report", "report", "Summarize a similarity matrix."),
    # plot
    "plot_contigs": Tool("plot_contigs", "plot", "Draw contigs and their annotations as an interactive HTML page."),
    "color_table_to_legend": Tool("color_table_to_legend", "plot", "Draw an SVG legend from a color table."),
    # database
    "db_download": Tool("db_download", "database", "Download sequence databases in formats Domainator reads."),
    "format_db": Tool("format_db", "database", "Shard, compress, and index a database for domain_search."),
}
