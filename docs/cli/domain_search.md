[CLI reference](README.md)
# domainator domain_search

Search a large database with a few references and extract the hits.

```
Search for matches to protein or nucleic-acid references.

For searching large databases of sequences (genbank files or fasta files) with small numbers (< ~100) of queries.
Protein references are searched against CDSs or protein records.
Nucleic-acid references are searched against whole nucleotide contigs.
Supported references include protein HMMs, protein FASTA, nucleotide FASTA, nucleotide HMMs, and infernal CM files.
Protein inputs cannot be searched with nucleotide references.

Returns the contigs with hits, possibly truncated to the neighborhood around the hits.

When --max_hits is specified, it applies to the entire search, not individual queries. Thus if you supply 5 query profiles and --max_hits 100, 
at most 100 hits will be returned (not 500). 

When max_hits or pad are specified, hits will be sorted by best domain score and written after the entire search completes.
If neither is set, then hits will be written on the fly and not sorted.
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | the genbank or fasta files to search. Genbank files can be nucleotide (with CDS annotations) or peptide. Fasta files can be protein or nucleotide; use --alphabet to specify which. |
| `--alphabet` | {nucleotide,protein} | protein | Whether the sequences in FASTA inputs are protein or nucleotide. GenBank records declare their own molecule type. |
| `-r`, `--reference` | FILE ... | required | reference files to search with. Supported types are protein HMMs, protein FASTA, nucleotide FASTA, nucleotide HMMs, and infernal CM files. |
| `-o`, `--genbank_out` | FILE |  | Write the hits to this GenBank file. If the filename ends in .gz/.bgz, the output is written as BGZF (block-gzip): it is gunzip-readable and ~3-4x smaller, and can be used directly as a random-access, partitionable compressed input to a later domainator step. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the hits to this FASTA file. Use - for standard output. |
| `-Z` | INTEGER | 1000 | Passed as the -Z parameter to hmmsearch: Assert that the total number of targets in your searches is <x>, for the purposes of per-sequence E-value calculations, rather than the actual number of targets seen.  Set to 0 to use the number of target sequences in the input, counted before any taxonomy filter, so filtering does not change E-values. Supplying -Z is recommended for most use cases of domainate, as it will speed up the search and make comparisons between searches more meaningful. |
| `-e`, `--evalue` | NUMBER | 0.001 | threshold E value for the domain hit. Must be <=10. |
| `--min_evalue` | NUMBER | 0.0 | hits with E value lower than this will be filtered out. Not frequently used. Use only if you want to eliminate close matches. Must be >=0. [default 0] |
| `--max_mode` |  |  | Run hmmsearch, phmmer, or nhmmer in maximum sensitivity mode, which is much slower, but more sensitive. |
| `--max_hits` | INTEGER |  | the maximum number of CDSs returned by the search, keeping the best-scoring ones by the bitscore of their best scoring profile. If not supplied, return all hits passing the evalue threshold. |
| `--max_hits_per_contig` | INTEGER |  | the maximum number of independent hits per contig. For nucleotide contig-level hits, each non-overlapping hit region is treated as an independent hit. [default: unlimited] |
| `--include_taxids` | INTEGER ... |  | Space separated list of taxids to include |
| `--exclude_taxids` | INTEGER ... |  | Space separated list of taxids to exclude |
| `--taxonomy_expr` | STRING |  | A boolean expression over taxids using operators & (AND), \| (OR), ~ (NOT), and parentheses, e.g. "2 & ~1224" (within Bacteria but not Proteobacteria). A taxid is true for a contig when it is in the contig's lineage. Mutually exclusive with --include_taxids/--exclude_taxids. |
| `--ncbi_taxonomy_path` | DIR | /tmp/ncbi_taxonomy_sean | Path to NCBI taxonomy database directory. Will be created and downloaded if it does not exist. |
| `--taxonomy_update` |  |  | If taxonomy database exists, check it against the version on the ncbi server and update if there is a newer version. |
| `--partial` | {include,exclude,only} | include | Which CDSs/proteins to search, by fragment status. include: all. exclude: skip fragments. only: only fragments. A CDS is a fragment if its location has '<' or '>', or it has a /partial qualifier. A protein record is a fragment if it has a UniProt ' (Fragment)' or MGnify 'FL=0' header, a 'Flags: Fragment' definition or NON_TER features (UniProt text entries), or a CDS/Protein feature with '<' or '>' spanning it. |
| `--decoys` | STRING ... |  | Names of decoy domains. A decoy domain is a domain that is not expected to be found in the input sequences. If the best hit for a CDS/protein is a domain from the decoys list, it will be not be returned as a hit. |
| `--decoys_file` | FILE |  | text file containing names of decoy domains. |
| `--cpu` | INTEGER | 0 | the number of cores of the cpu which are used at a time to run the search (0: one thread per physical core) |
| `--max_overlap` | NUMBER | 1 | the maximum fractional overlap between domains to be included in the annotated genbank. If >= 1, then no overlap filtering will be done. |
| `--add_annotations` |  |  | When activated, new domainator annotations will be added to the file, not just the domain_search annotations. Useful if you want to see what the non-best hits score. |
| `--deduplicate`, `--no_deduplicate` |  |  | if the same region is extracted for multiple hits, then keep only one copy of it (compared by record ID). |
| `--deduplicate_by_sequence` |  |  | Drop hits whose output sequence is identical to an earlier hit (SHA-256 hash, as in domainator deduplicate_genbank --algorithm hash). Useful for highly redundant databases. Independent of --deduplicate, which only drops repeated record IDs. |
| `--both_strands` |  |  | With --deduplicate_by_sequence, also treat complement and reverse-complement nucleotide sequences as duplicates. |
| `--max_region_overlap` | NUMBER | 1.0 | the maximum fractional of overlap between any two output regions. If >= 1, then no overlap filtering will be done. Regions are output in a greedy fashion based on CDS start site. New regions are output if less than this fraction of them overlaps with any previously output region. |
| `--translate` |  |  | by default, nucleotide inputs return nucleotide regions. When --translate is set, CDS hits are translated before writing. Whole-contig nucleic-acid hits are not translated. Note that --translate is incompatble with neighborhood extraction. |
| `--kb_range` | NUMBER |  | How many kb upstream and downstream of the selected CDS to extract. Starting from the ends of the CDS. |
| `--kb_up` | NUMBER |  | How many kb upstream of the selected CDS to extract. Starting from the end of the CDS. |
| `--kb_down` | NUMBER |  | How many kb downstream of the selected CDS to extract. Starting from the end of the CDS. |
| `--cds_range` | INTEGER |  | How many CDSs upstream and downstream of the selected CDS to extract. Example, 1 would mean to return contigs including the selected CDSs and the immediate upstream and downstream CDSs. |
| `--cds_up` | INTEGER |  | How many CDSs upstream of the selected CDS to extract. Example, 1 would mean to return contigs including the selected CDSs and the immediate upstream CDS. |
| `--cds_down` | INTEGER |  | How many CDSs downstream of the selected CDS to extract. Example, 1 would mean to return contigs including the selected CDSs and the immediate downstream CDS. |
| `--whole_contig` |  |  | If set, then return the entire matching contigs. Circular contigs stay circular, and are rotated so that the hit is in the center. |
| `--strand` | {f,r} |  | Only extract regions around CDSs on the specified strand. |
| `--pad` |  |  | If set, then ends of sequences will be padded so that the focus CDS aligns. See --anchor. |
| `--anchor` | {start,end,center} | start | Where to align the focus CDSs when --pad is set, relative to the focus CDS in its own orientation. start: align the first base of the focus CDSs. end: align their last base. center: align their midpoints. |
| `--keep_direction` |  |  | by default extracted regions will be flipped so that the focus cds is on the forward strand. Setting this option will keep the focus cds on whatever strand it started on. |
| `--batch_size` | INTEGER | 10000 | Approximately how many target sequences to search at one time in a batch. This also controls partition size: each parallel task processes about this many CDSs/sequences. Larger values mean fewer, bigger partitions (lower per-task overhead, but coarser load balancing); smaller values give finer parallelism at the cost of more tasks. The input is scanned once to find partition boundaries regardless of this value. |
| `--gene_call` | {all,unannotated} |  | When activated, new CDS annotations will be added with Prodigal in Metagenomic mode. If 'all', then any existing CDS annotations will be deleted and all contigs will be re-annotated. If 'unannotated', then only contigs without CDS annotations will be annotated. -Z 0 can't be used with gene calling, because the number of target sequences isn't known until the genes are called. |
| `--max_output_gb` | NUMBER | 25.0 | Fail if a projected output file would exceed this many GB. Set to 0 to disable the guardrail. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
