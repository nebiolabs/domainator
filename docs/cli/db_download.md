[CLI reference](README.md)
# domainator db_download

Download sequence databases in formats Domainator reads.

```
Download sequence databases from the internet in formats appropriate for Domainator

Including the following databases:
- Uniprot (includes Trembl and SwissProt) (downloads fasta by default, use _gb suffix to download genbank file with addtional annotations)
- Trembl (downloads fasta by default, use _gb suffix to download genbank file with addtional annotations)
- SwissProt (downloads fasta by default, use _gb suffix to download genbank file with addtional annotations)
- ncbi_complete_genome_proks: NCBI genomes where assembly_level is "Complete Genome". (excluding eukaryotes)
- ncbi_complete_genome_nonredundant_proks: non-redundant NCBI genomes where assembly_level is "Complete Genome". (excluding eukaryotes)
- ncbi_representative_proks : representative and reference genomes from the NCBI database (excluding eukaryotes)
- ncbi_representative_all : representative and reference genomes from the NCBI database (all taxa, including eukaryotes)
- ncbi_nonredundant_proks : non-redundant genomes from the NCBI database (excluding eukaryotes)
- ncbi_nonredundant_all : non-redundant genomes from the NCBI database (all taxa, including eukaryotes)
- ncbi_all : all genomes from the NCBI database including duplicates of the same species (all taxa, including eukaryotes). Can produce very large files, up to 10 TB or more, so use with caution.


Key taxa:
    81490 : prokaryotic environmental samples
    408169 : metagenomes
    2 : bacteria
    2157 : archaea
    10239 : viruses
    
    2759 : eukaryota

WARNING: Downloading large databases can take a long time and use a lot of disk space!
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `--db` | {ncbi_all,ncbi_complete_genome_nonredundant_proks,ncbi_complete_genome_proks,ncbi_nonredundant_all,ncbi_nonredundant_proks,ncbi_representative_all,ncbi_representative_proks,swissprot,swissprot_gb,trembl,trembl_gb,uniprot,uniprot_gb} | required | Database to download. swissprot, trembl and uniprot are downloaded as FASTA (write them with --fasta_out); the _gb UniProt databases and the ncbi databases are downloaded as GenBank (write them with -o/--genbank_out). |
| `-o`, `--genbank_out` | FILE |  | Write the downloaded GenBank records to this file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the downloaded FASTA records to this file. Use - for standard output. |
| `--partial` | {include,exclude,only} | include | Which UniProt records to download, by fragment status. include: all. exclude: skip fragments. only: only fragments. A CDS is a fragment if its location has '<' or '>', or it has a /partial qualifier. A protein record is a fragment if it has a UniProt ' (Fragment)' or MGnify 'FL=0' header, a 'Flags: Fragment' definition or NON_TER features (UniProt text entries), or a CDS/Protein feature with '<' or '>' spanning it. Only for the UniProt databases. |
| `--include_taxids` | INTEGER ... |  | Space separated list of taxids to include |
| `--exclude_taxids` | INTEGER ... |  | Space separated list of taxids to exclude |
| `--taxonomy_expr` | STRING |  | A boolean expression over taxids using operators & (AND), \| (OR), ~ (NOT), and parentheses, e.g. "2 & ~1224" (within Bacteria but not Proteobacteria). A taxid is true for a record when it is in the record's lineage. Mutually exclusive with --include_taxids/--exclude_taxids. For prokaryote-only (_proks) databases, Eukaryota are excluded in addition to this expression. |
| `--ncbi_taxonomy_path` | DIR | /tmp/ncbi_taxonomy_sean | Path to NCBI taxonomy database directory. Will be created and downloaded if it does not exist. |
| `--taxonomy_update`, `--no_taxonomy_update` |  |  | if the taxonomy database exists, check it against the version on the ncbi server and update it if there is a newer version. |
| `--success_log_out` | FILE |  | Path to log file for successfully written records. If not specified, successfully written records will not be logged. |
| `--exclude_accessions` | STRING ... |  | Exclude genbank records with these accessions. Should be versionless (i.e. without the '.' and number after the '.') |
| `--exclude_accessions_file` | FILE |  | text file containing accessions genbank records to exclude. Should be versionless (i.e. without the '.' and number after the '.') |
| `--ncbi_summary` | FILE |  | Path to NCBI summary file. Normally it will be downloaded from the Internet on every execution, but you may want to pre-download it, or make custom subsets, the url is: https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt. |
| `--gene_call` | {all,unannotated} |  | When activated, for nucleotide databases, such as Genbank, new CDS annotations will be added with Prodigal in Metagenomic mode. If 'all', then any existing CDS annotations will be deleted and all contigs will be re-annotated. If 'unannotated', then only contigs without CDS annotations will be annotated. |
| `--max_records` | INTEGER |  | Stop after downloading this many records. For UniProt, this is the number of proteins. For genbank, the number of genomes, including multi-contig genomes. This is mainly for testing. |
| `--append` |  |  | Append to existing file |
| `--before` | STRING |  | For ncbi databases, only download records before this date. Format is YYYY-MM-DD. |
| `--after` | STRING |  | For ncbi databases, only download records after this date. Format is YYYY-MM-DD. |
| `--cpu` | INTEGER | 0 | the number of cores of the cpu which are used at a time to download ncbi records. 0: one per physical core. |
| `--download_backend` | {direct,datasets} | datasets | Network backend for ncbi genome downloads. 'datasets' uses the NCBI Datasets CLI dehydrate/rehydrate workflow (recommended for large jobs; NCBI's own bulk tool handles throttling). 'direct' downloads .gbff.gz files over HTTPS directly. |
| `--download_workers` | INTEGER | 3 | For the 'direct' backend: number of simultaneous HTTPS connections to the NCBI genomes server. NCBI rate-limits by concurrency, so keep this low (~3). 0 = use --cpu. With --gene_call set, prodigal runs inside these workers. |
| `--email` | STRING |  | Contact email included in the User-Agent header sent to NCBI (politeness; not an API key). Falls back to the NCBI_EMAIL environment variable. |
| `--ncbi_api_key` | STRING |  | NCBI API key passed to the 'datasets' backend (does not affect the raw FTP server / 'direct' backend). Falls back to the NCBI_API_KEY environment variable. |
| `--datasets_workdir` | DIR |  | For the 'datasets' backend: working directory for the dehydrated package and rehydrated files. Required when --download_backend is 'datasets'. |
| `--datasets_include` | STRING | gbff | For the 'datasets' backend: value passed to `datasets ... --include`. Only the GenBank flat file (gbff) is used downstream, so the genomic FASTA is not requested by default. \[default: gbff\] |
| `--datasets_max_workers` | INTEGER | 10 | For the 'datasets' backend: value passed to `datasets rehydrate --max-workers` (allowed range 1-30). 10 is the NCBI Datasets default and matches the 10 req/s API-key limit. \[default: 10\] |
| `--datasets_path` | PROGRAM |  | Path to the `datasets` (NCBI Datasets CLI) binary. If not given, it is looked up on PATH. |
| `--datasets_gzip` |  |  | For the 'datasets' backend: pass `--gzip` to `datasets rehydrate` so the rehydrated genomic.gbff files are stored gzip-compressed (much smaller scratch dir; read transparently). Recommended for very large downloads. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
