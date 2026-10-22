[CLI reference](README.md)
# domainator domainate

Annotate contigs or proteins with hits to HMM, sequence, or CM references.

```
Annotate sequence files with protein or nucleic-acid references.

Protein references are searched against CDSs or protein records with hmmsearch or phmmer.
Nucleic-acid references are searched against whole nucleotide contigs with nhmmer or infernal cmsearch.
Hits passing the threshold are added to the output GenBank file as new Domainator features.

Supported references include protein HMMs, protein FASTA, nucleotide FASTA, nucleotide HMMs, and infernal CM files.
Protein inputs cannot be searched with nucleotide references.

NOTES: genbank files that contain translation annotations in their CDS features will run much faster than those without translation annotations.
    domainator domainate stores the entire reference set in memory, so it is not suitable for cases where the reference
    set is larger than your system memory (this may change in future versions).
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | the genbank or fasta files to annotate. Genbank files can be nucleotide (with CDS annotations) or peptide. Fasta files can be protein or nucleotide; use --alphabet to specify which. |
| `--alphabet` | {nucleotide,protein} | protein | Whether the sequences in FASTA inputs are protein or nucleotide. GenBank records declare their own molecule type. |
| `-r`, `--reference` | FILE ... |  | reference files to search against the input. Supported types are protein HMMs, protein FASTA, nucleotide FASTA, nucleotide HMMs, and infernal CM files. |
| `--foldseek` | FILE ... |  | Foldseek database files. |
| `--esm2_3di_weights` | FILE |  | checkpoint file for esm2_3Di model. |
| `--esm2_3di_device` | STRING | cuda:0 | device to use for esm2_3Di model. |
| `--foldseek_device` | STRING | cpu | device to use for the foldseek search. 'cpu' for a CPU-only search, or a CUDA device string ('cuda', 'cuda:0', 'cuda:1', ...) to enable GPU-accelerated search. GPU search requires a recent foldseek and a GPU-compatible (padded) target database (see 'foldseek makepaddedseqdb'). |
| `--foldseek_path` | PROGRAM |  | path to the foldseek executable, if it is not on the PATH. |
| `-o`, `--genbank_out` | FILE |  | Write the records to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the records to this FASTA file. Use - for standard output. |
| `-Z` | INTEGER | 1000 | Passed as the -Z parameter to hmmsearch: Assert that the total number of targets in your searches is <x>, for the purposes of per-sequence E-value calculations, rather than the actual number of targets seen.  Set to 0 to use the number of target sequences in the input, counted before any taxonomy filter, so filtering does not change E-values. Supplying -Z is recommended for most use cases of domainate, as it will speed up the search and make comparisons between searches more meaningful. |
| `-e`, `--evalue` | NUMBER | 0.001 | threshold E value for the domain hit. Must be <=10. |
| `--min_evalue` | NUMBER | 0.0 | hits with E value LOWER (as in BETTER) than this will be filtered out. NOT FREQUENTLY USED. Use only if you want to eliminate close matches. Must be >=0. \[default 0\] |
| `--max_hits_per_cds` | INTEGER |  | the maximum number of domains to be reported per CDS. If not specified, then all hits are reported. |
| `--max_mode` |  |  | Run hmmsearch, phmmer, or nhmmer in maximum sensitivity mode, which is much slower, but more sensitive. |
| `--include_taxids` | INTEGER ... |  | Space separated list of taxids to include. Contigs with taxonomy not in this list will be skipped. |
| `--exclude_taxids` | INTEGER ... |  | Space separated list of taxids to exclude. Contigs with taxonomy in this list will be skipped. |
| `--taxonomy_expr` | STRING |  | A boolean expression over taxids using operators & (AND), \| (OR), ~ (NOT), and parentheses, e.g. "2 & ~1224" (within Bacteria but not Proteobacteria). A taxid is true for a contig when it is in the contig's lineage. Mutually exclusive with --include_taxids/--exclude_taxids. |
| `--partial` | {include,exclude,only} | include | Which CDSs/proteins to search, by fragment status. include: all. exclude: skip fragments. only: only fragments. A CDS is a fragment if its location has '<' or '>', or it has a /partial qualifier. A protein record is a fragment if it has a UniProt ' (Fragment)' or MGnify 'FL=0' header, a 'Flags: Fragment' definition or NON_TER features (UniProt text entries), or a CDS/Protein feature with '<' or '>' spanning it. |
| `--ncbi_taxonomy_path` | DIR | /tmp/ncbi_taxonomy_sean | Path to NCBI taxonomy database directory. Will be created and downloaded if it does not exist. |
| `--taxonomy_update` |  |  | If taxonomy database exists, check it against the version on the ncbi server and update if there is a newer version. |
| `--cpu` | INTEGER | 0 | the number of cores of the cpu which are used at a time to run the search (0: one thread per physical core) |
| `--max_overlap` | NUMBER | 1 | the maximum fractional of overlap between domains to be included in the annotated genbank. If >= 1, then no overlap filtering will be done. |
| `--overlap_by_db` |  |  | If activated, then overlap filtering will be done by database, rather than all together. |
| `--hits_only` |  |  | when activated, the ouptut will only have contigs with at least one domain annotation. In many cases, domainator domain_search will be faster and more appropriate. |
| `--no_annotations` |  |  | when activated, new annotations will not be added to the file. This is useful if you are just trying to extract a set of contigs for annotation in a later step. |
| `--batch_size` | INTEGER | 10000 | How many target sequences to search at one time in a batch. Increasing this number might improve speeds at the cost of memory. |
| `--gene_call` | {all,unannotated} |  | When activated, new CDS annotations will be added with Prodigal in Metagenomic mode. If 'all', then any existing CDS annotations will be deleted and all contigs will be re-annotated. If 'unannotated', then only contigs without CDS annotations will be annotated. If you supply a nucleotide fasta file, be sure to also specify --alphabet nucleotide. -Z 0 can't be used with gene calling, because the number of target sequences isn't known until the genes are called. |
| `--max_output_gb` | NUMBER | 25.0 | Fail if a projected output file would exceed this many GB. Set to 0 to disable the guardrail. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
