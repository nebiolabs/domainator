[index](README.md)
# Structure search and annotation

Domainator can use protein structures as input, either annotating them with reference
structures or searching a large structure database with a few reference structures. Both
produce ordinary Domainator genbank files, so the rest of the suite — `extract_domains.py`,
`plot_contigs.py`, `enum_report.py`, `summary_report.py`, `color_genbank.py`,
`build_ssn.py`, `seq_dist.py`, `select_by_contig.py` — works on the results unchanged.

There are three tools, mirroring the sequence-based ones:

|  | sequence input | structure input |
| ---- | ---- | ---- |
| **annotate** (all records, `Domainator` features) | `domainate.py` | `structure_domainate.py` |
| **search** (hit records, `Domain_Search` features) | `domain_search.py` | `structure_search.py` |
| **convert** | — | `structure_to_genbank.py` |

## Requirements

The structure tools shell out to [foldseek](https://github.com/steineggerlab/foldseek),
which is installed with the conda environment. Nothing else is needed — in particular these
tools do **not** need PyTorch, `esmologs`, or the ESM-2 3B checkpoint. That stack is only
needed for [annotating sequences that have no structure](esm_3b_foldseek.md).

Check that foldseek is available with `foldseek version`, or pass an explicit path with
`--foldseek_path`.

## Accepted inputs

Both `-i` and `-r` accept, in any of the three tools:

* individual structure files — `.pdb`, `.cif`, `.mmcif`, `.ent`, `.bcif`, optionally gzipped,
* a directory, searched recursively for structure files,
* a glob such as `structures/*.pdb.gz` (expanded internally, so it also works from a
  `--config` file, where the shell would not expand it),
* the prefix of a prebuilt foldseek database.

One record is written per structure **chain**. Chain splitting, structure parsing, and entry
naming are all done by foldseek, so records are named exactly as they appear in a hand-run
foldseek search: `<file basename>_<chain>`, for example `1ubq_A`.

## Annotating structures

```bash
structure_domainate.py -i my_structures/ -r reference_structures/ -o annotated.gb
plot_contigs.py -i annotated.gb -o annotated.html
```

Every input chain is written, whether or not it had a hit. Add `--hits_only` to keep just
the chains that matched.

Because both databases are built from real structures they carry C-alpha coordinates, so the
full range of foldseek alignment modes is available:

```bash
structure_domainate.py -i my_structures/ -r reference_structures/ \
    --alignment_type 1 --metrics tmscore lddt rmsd -o annotated.gb
```

`--alignment_type` is `0` for 3Di only, `1` for TM-align, and `2` (the default) for 3Di+AA.
`--metrics` adds `tmscore`, `lddt`, `rmsd` and `prob` qualifiers to the annotations.
`tmscore`, `lddt` and `rmsd` need coordinates in both databases; asking for them against a
database built from sequences rather than structures is an error rather than a silent
degradation.

## Searching a large structure database

```bash
structure_search.py -i /path/to/afdb_foldseek_db -r query.pdb \
    --max_hits 1000 --max_seqs 1000000 --device cuda:0 --cpu 6 -o hits.gb
enum_report.py -i hits.gb --domains -o hits.tsv
```

Only records that matched are written, each carrying a `Domain_Search` feature marking the
region homologous to its best-matching reference. `--add_annotations` also writes a
`Domainator` feature per hit.

Hit records are built from the sequences foldseek reports for each hit, so the searched
database is never read end to end and its size does not bound memory usage. Hits are grouped
by database record using an external sort, so memory scales with the hits for one record
rather than with the whole result set.

Two knobs matter at scale:

* `--max_seqs` bounds how many hits foldseek keeps **per reference**. The default of 1000 is
  low for a database of hundreds of millions of structures. A warning is printed whenever a
  reference saturates it, so truncation is never silent.
* `--max_hits` bounds how many **records** are returned, keeping the best-scoring ones. It
  applies to the whole search, not per reference, and makes the output sorted best-first and
  written after the search completes.

## Reusing a database

Building a structure database is the expensive part of a run, so it can be persisted and
reused:

```bash
# build the reference database once
structure_to_genbank.py -i reference_structures/ --keep_db refdb -o refs.gb

# reuse it, and keep the input database too
structure_domainate.py -i my_structures/ -r refdb --keep_db mydb -o annotated.gb
structure_search.py    -i mydb -r refdb -o hits.gb
```

A database written by `--keep_db` carries C-alpha coordinates, so it supports every alignment
mode and metric. Databases built from amino acid and 3Di sequences instead — including the
ones described in the [ESM-2 3B 3Di documentation](esm_3b_foldseek.md) — do not, and are
limited to `--alignment_type 0` and `2`.

## Using structures with the sequence-based tools

`structure_to_genbank.py` converts structures without searching, which makes them ordinary
Domainator inputs:

```bash
structure_to_genbank.py -i my_structures/ -o chains.gb

domainate.py -i chains.gb -r pfam.hmm -o annotated.gb        # annotate with HMM profiles
domain_search.py -i chains.gb -r query.hmm -o hits.gb        # search them by sequence
seq_dist.py -i chains.gb -r chains.gb --sparse distances.hdf5
deduplicate_genbank.py -i chains.gb -o nonredundant.gb
```

All three structure tools write to stdout when `-o` is omitted or given as `-`, so they can
be piped into any tool that reads genbank from stdin:

```bash
structure_to_genbank.py -i my_structures/ -o - | color_genbank.py --color_domains -o colored.gb
```

Note that `domainate.py` and `domain_search.py` read from files rather than stdin, so write
an intermediate file for those (as above).

`--store_3di` additionally records each chain's 3Di structural-alphabet string in a
`/threedi` qualifier, which is useful for rebuilding a foldseek database from a genbank file.
Note that it is one opaque string per record, so it does not survive tools that slice
records.

## Choosing an alignment direction, and what `--evalue` means

References are always aligned as the **query** and the input as the **target**. This is what
foldseek is optimized for, and it means `--evalue` has the same meaning it has everywhere
else in Domainator: it is computed against the database being searched or annotated. Running
it the other way would make `--evalue` depend on how many reference structures you happened
to supply.

One consequence worth knowing: E-values shift when the input database changes size, so the
same reference against the same structure will score differently in a 9-record and a
9-million-record database. This is the same behavior as `hmmsearch`, where `-Z` exists to
normalize it; foldseek has no equivalent knob.

Annotation coordinates are positions on the **input** record. The matching positions on the
reference are recorded in the `/rstart`, `/rend` and `/rlen` qualifiers. See
[File Formats](file_formats.md) for the full qualifier list.

`--min_evalue` excludes hits at or below a given E-value, which is the usual way to drop a
reference's match to its own copy inside the database. foldseek has no equivalent option, so
Domainator applies this bound itself.

## Cross-checking against foldseek directly

`--hits_tsv` writes the raw foldseek hit table alongside the genbank output, which makes the
results directly comparable to a hand-run search:

```bash
foldseek easy-search refs/*.pdb $db ref.tsv tmp --format-mode 4 -e 0.001 \
    --format-output query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen
structure_search.py -i $db -r refs/ -e 0.001 --hits_tsv ours.tsv -o hits.gb
```

## Passing options through to foldseek

`--aligner_arg` appends arguments verbatim to the aligner's search command line, for options
Domainator does not wrap:

```bash
structure_search.py -i $db -r refs/ --aligner_arg --prefilter-mode 2 -o hits.gb
```

## Limitations

* **Residue numbering.** Coordinates are indices into the observed residue sequence foldseek
  extracts from the structure, not original PDB residue numbers. A chain with a disordered
  loop has a sequence shorter than its residue-number span, so hit positions cannot be
  mapped back to the structure's own numbering.
* **No coordinates in genbank.** Structure coordinates are deliberately not stored in the
  genbank output — they would inflate records several-fold, exceed what standard genbank
  parsers accept on one line, and would be silently wrong after any tool that slices
  records. Persist the aligner database with `--keep_db` instead.
* **One backend.** foldseek is currently the only backend. `--algorithm` and `--aligner_arg`
  exist so that others can be added without changing the interface.
