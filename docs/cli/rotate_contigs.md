[CLI reference](README.md)
# domainator rotate_contigs

Set the origin of circular contigs relative to a domain.

```
Rotate circular contigs so that their origin is set relative to a domain annotation.

For each contig, one domain is chosen from new annotations found by searching the references (-r) against the contig,
and/or from pre-existing Domainator (--domains) or domain_search (--annotation_type domain_search) annotations.
If there are several candidates, the one with the highest score is chosen (ties broken by lowest E-value, then by position).
The contig is then rotated so that the origin is at the start, end, or center of the chosen domain (--anchor),
and reverse complemented if necessary to put the chosen domain on the forward strand (unless --keep_direction).

Positions are relative to the chosen domain, in its own orientation:
    start:  the first base of the domain becomes the first base of the contig.
    end:    the last base of the domain becomes the last base of the contig.
    center: the origin is placed in the middle of the domain.

Inputs must be nucleotide sequences. Only circular contigs are rotated. FASTA records and linear GenBank records are
skipped with a warning, unless --circularize is set, in which case they are treated as circular.
Contigs without a qualifying domain are skipped with a warning. Output contigs are always marked as circular.

Protein references are searched against CDS translations, so contigs without CDS annotations (such as all FASTA inputs)
need --gene_call to be found by protein references. Nucleotide references (nucleotide HMMs, nucleotide FASTA, and
infernal CM files) are searched against the whole contig.
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | Names of input genbank or fasta files. Must be nucleotide sequences. Use - to read standard input. |
| `-o`, `--genbank_out` | FILE |  | Write the records to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the records to this FASTA file. Use - for standard output. |
| `-r`, `--reference` | FILE ... |  | reference files to search against the input, to find new domains to rotate on. Supported types are protein HMMs, protein FASTA, nucleotide FASTA, nucleotide HMMs, and infernal CM files. |
| `-Z` | INTEGER | 1000 | Passed as the -Z parameter to hmmsearch: Assert that the total number of targets in your searches is <x>, for the purposes of per-sequence E-value calculations, rather than the actual number of targets seen. Set to 0 to use the number of actual target sequences in the input files (not available when reading from stdin). |
| `-e`, `--evalue` | NUMBER | 0.001 | threshold E value for domains to rotate on. Both new and pre-existing annotations must pass this threshold. Must be <=10. |
| `--domains` | STRING ... |  | Rotate on pre-existing domainator annotations with these names. |
| `--domains_file` | FILE |  | text file containing names of pre-existing domainator annotations to rotate on. |
| `--databases` | STRING ... |  | Consider only pre-existing annotations from these databases. default: all databases. |
| `--annotation_type` | {domainator,domain_search} | domainator | Which pre-existing annotations to rotate on. domainator: domain annotations (from domainate). domain_search: best-hit annotations (from domain_search). |
| `--anchor` | {start,end,center} | start | where to put the origin, relative to the chosen domain in its own orientation (not relative to the original orientation of the contig). start: the first base of the domain becomes the first base of the contig. end: the last base of the domain becomes the last base of the contig. center: the origin is placed in the middle of the domain. |
| `--keep_direction` |  |  | by default, contigs are reverse complemented if necessary to put the chosen domain on the forward strand. Setting this option will keep the original orientation of the contig. |
| `--circularize` |  |  | treat fasta records and linear genbank records as circular, before searching and rotating. By default, they are skipped. |
| `--keep_name` |  |  | by default, rotated contigs will be renamed to note the original coordinates of their first and last bases, and whether they were reverse complemented. Set this option to keep contigs with their original names. |
| `--keep_annotation` |  |  | by default, new annotations from the reference search are not written to the output. Setting this option will keep the new annotation of the domain that the contig was rotated on. Has no effect if the contig was rotated on a pre-existing annotation. |
| `--gene_call` | {all,unannotated} |  | When activated, new CDS annotations will be added with Prodigal in Metagenomic mode before the search. If 'all', then any existing CDS and domainator annotations will be deleted and all contigs will be re-annotated. If 'unannotated', then only contigs without CDS annotations will be annotated. Needed to search protein references against contigs without CDS annotations, such as fasta records. |
| `--cpu` | INTEGER | 0 | the number of cores of the cpu which are used at a time to run the search (0: one thread per physical core) |
| `--max_mode` |  |  | Run hmmsearch, phmmer, or nhmmer in maximum sensitivity mode, which is much slower, but more sensitive. |
| `--batch_size` | INTEGER | 1000 | How many contigs to search at one time in a batch. Increasing this number might improve speeds at the cost of memory. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
