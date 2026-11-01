# `db_download` fixtures

Local mirrors that let `test_db_download.py` run `db_download` offline, through the
same code paths as a real download.

- `uniprot/` stands in for UniProt's `knowledgebase/complete` release directory (read
  with `--uniprot_dir`). Its files hold real UniProt entries, copied from the
  2026-11 release (UniProt is CC BY 4.0): SwissProt C0HLM9 and P80438 (fragments),
  P0C6K8 (a virus, with host lines), Q01457 (from a plasmid) and P0A7Y4 (many
  features and references); TrEMBL A0A0U3GSY3, Q50831 (a fragment) and
  A0A0K2RVI7. The FASTA files hold the same entries.
- `taxdmp/` is a small NCBI taxonomy covering those entries' taxids and the
  assemblies'. Lineages are shortened to the domain.
- `assembly_summary.txt` (for `--ncbi_summary`) and `genomes/` stand in for NCBI's
  assembly summary and genomes server. The `ftp_path` column names directories under
  `genomes/`, relative to the summary, whose `*_genomic.gbff.gz` files are gzipped
  copies of other fixtures. GCA_000000088's is a truncated download,
  GCA_000000099 has no directory and GCA_000000101 has no `ftp_path`.
- `exclude_accessions.txt` is for `--exclude_accessions_file`.
