[CLI reference](README.md)
# domainator extract_domains

Extract annotated domains as their own records.

```
Extract domain-annotated regions

Takes a genbank file which has been annotated using the domainator program, looks for the specified domains, and extracts them from the contigs.

This tool is typically used with protein sequences, but can also be used with nucleic acid sequences.
Handling of nucleic acids is less stable and will may change in future updates, so use with caution.

can write a fasta or a genbank file
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | Genbank filenames. Use - to read standard input. |
| `-o`, `--genbank_out` | FILE |  | Write the records to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the records to this FASTA file. Use - for standard output. |
| `--contigs` | STRING ... |  | only extract domains from contigs with ids in this list. Additive with --contigs_file. default: all contigs. |
| `--contigs_file` | FILE |  | only extract domains from contigs with ids listed in this file (one per line). Additive with --contigs. default: all contigs. |
| `-e`, `--evalue` | NUMBER |  | the evalue cutoff (max) for domain annotations to be extracted |
| `--min_score` | NUMBER |  | the score cutoff (min) for domain annotations to be extracted |
| `--domains` | STRING ... |  | Include domains with these names |
| `--domains_file` | FILE |  | text file containing names of domains to include. |
| `--databases` | STRING ... |  | Domain databases to use for domain extraction. default: all databases. |
| `--annotation_type` | {domainator,domain_search} | domainator | Which annotations to extract. domainator: domain annotations (from domainate). domain_search: best-hit annotations (from domain_search). |
| `--partial` | {include,exclude,only} | include | Which domains to extract, by the CDS or protein they are on, by fragment status. include: all. exclude: skip fragments. only: only fragments. A CDS is a fragment if its location has '<' or '>', or it has a /partial qualifier. A protein record is a fragment if it has a UniProt ' (Fragment)' or MGnify 'FL=0' header, a 'Flags: Fragment' definition or NON_TER features (UniProt text entries), or a CDS/Protein feature with '<' or '>' spanning it. Contig-level nucleotide hits are never on fragments. |
| `--extend_range` | INTEGER |  | extract this many additional residues (or bases) on each side of the hit region. Sets the default of --extend_up and --extend_down. |
| `--extend_up` | INTEGER |  | extract this many additional residues (or bases) upstream of the hit region. Default: --extend_range, or 0. |
| `--extend_down` | INTEGER |  | extract this many additional residues (or bases) downstream of the hit region. Default: --extend_range, or 0. |
| `--combine` |  |  | if set, then when multiple target domains occur in the same contig, an envelope containing all of the hit domains will be extracted. By default, each hit domain will be extracted separately. This option will raise an exception on circular contigs, because it would be ambiguous. Cannot be used with --splice. |
| `--splice` |  |  | if set, then when multiple target domains occur in the same contig, the intermediate sequence will be deleted and the extracted domains from the contig glued together, in the order they are read on the strand of the domain with the lowest start coordinate. By default, each hit domain will be extracted separately. Cannot be used with --combine. |
| `--keep_direction` |  |  | by default extracted regions will be flipped so that the domain (or the domain with the lowest start coordinate when --combine or --splice is active) is on the forward strand. Set this option to not reverse complement anything. In the rare case of a split domain with parts on different strands, the strand of the first part will be considered. |
| `--keep_name` |  |  | by default, extracted regions will be renamed to note their original coordinates and whether they were reverse complemented. Set this option to keep contigs with their original names. WARNING: if you aren't careful, this can easily result in duplicate names. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
