[index](README.md)
# Domainator file formats

## Genbank

Genbank files are the sequence format most used by Domainator. Genbank files are convenient because they can store information on both sequence and sequence positional annotations.

### Features
Domainator produces and understands some feature types particular to it.

`Domainator` features are domain annotations added to genbank files by `domainate.py` or `domain_search.py`. The `cds_id` qualifier is particularly critical as it links the `Domainator` feature to a `CDS` feature. Domainator adds a `cds_id` qualifier to any `CDS` feature lacking one. The `cds_id` is expected to be unique within a contig. The general format for an automatically generated `cds_id` is `[start coordinate]_[strand]_[end coordinate]`. Multi-part CDSs, for example CDSs containing introns or crossing the origin are handled in a somewhat more complicated way, as special cases, with the goal being to maintain a one-cds_id to one_cds relation. For nucleotide contig-level annotations (produced by nhmmer or infernal searches), `cds_id` is set to `"."` because the hit is not associated with a particular CDS.

`rstart`, `rend` are the start and end coordinates of the local alignment in the "reference" (i.e. profile). `rlen` is the total length of the reference profile, so you can get an idea of how the local alignment is oriented with respect to the entire reference profile/sequence.

```javascript
     Domainator      560..1165
                     /program="hmmsearch"
                     /database="pdonr_hmms"
                     /description="Chloramphenicol acetyltransferase"
                     /evalue="4.0e-100"
                     /score="329.7"
                     /name="CAT"
                     /identity="55.9"
                     /cds_id="2265_-1_1606"
                     /rstart="1"
                     /rend="100"
                     /rlen="100"
```


`Domain_Search` features are added by `domain_search.py`. `Domain_Search` annotations are related to `Domainator` annotations in that they have exactly the same qualifiers. What makes them distinct is that they are cleared from sequences used as input to `domain_search.py`, and there are exactly one per contig on sequences returned by `domain_search.py`. Other programs, such as `extract_domains.py`, `select_by_cds.py`, and others also handle `Domain_Search` annotations distinctly from `Domainator` annotations, typically by using the `--search_hits` command line argument.

### Taxonomy

Some Domainator scripts recognize Taxids from the [NCBI taxonomy database](https://www.ncbi.nlm.nih.gov/taxonomy). Domainator understands several ways of marking taxonomy. 

If a contig contains ` OX=[integer]` in its description, the integer will be understood as a taxid. This is common in Fasta files From UniProt, including SwissProt.

If that is not present, the `source` features will be examined in order from longest to shortest until a source annotation with `db_xref="taxon:[integer]"` is found. For example:

```javascript
     source          1..256
                     /db_xref="taxon:654924"
```

If no taxid is noted in either place, the Taxid will be assigned [32644](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=32644) for "unidentified". This taxid will be used internally by Domainator, but will not be written to any output files.


### Differences from BioPython in internal storing of Genbank records

`domainator.utils.parse_seqfiles` always assigns a `molecule_type` attribute to sequence records. For genbank files, it also swaps the `record.name` field and `record.id` field before returning the record. The `write_genbank` function swaps them back before writing. There are also some bugfixes in the parser relative to the Biopython parser to handle weird edge-case genbank files we've run into.

```python 
    id = record.id
    record.id = record.name
    record.name = id
```

## hdf5

Domainator uses hdf5 files to store data matrices in both dense and sparse formats. Code for reading, writing, and using data matrices is found in the `data_matrix.py` file (which is not executable, but imported by various other scripts).

hdf5 output is generally specified using the `--dense [filename]` or `--sparse [filename]` command line arguments, for example by `seq_dist.py`. For matrix input, the data format is automatically detected.

Domainator data matrix hdf5 files have the following attributes:
```python
ARRAY_TYPE: {DENSE, SPARSE_CSR} # indicates whether the hdf5 file describes a dense or sparse matrix. Sparse matrices are stored in Compressed Sparse Row format.
SYMMETRIC_LABELS: {True, False} # indicates whether the x and y axis labels are the same
MATRIX_FILE_VERSION: str # will be incremented on any backwards-compatibility-breaking changes to the matrix format.
(optional) DATA_TYPE: str # describes the type of data in the matrix (e.g. 'score', 'norm_score', 'row_norm_score', 'score_dist', 'bool', 'efi_score')
```

Domainator data matrix hdf5 files have the following datasets:
```python
ROW_LABELS: list of str
COL_LABELS: list of str # not used if SYMMETRIC_LABELS is True
DENSE_DATA: 2d array of float # used for dense data
SPARSE_VALUES: list of float # used for values for CSR sparse data
SPARSE_CSR_INDICES: list of int # column indices for CSR sparse data
SPARSE_CSR_INDPTR: list of int # row pointer offsets for CSR sparse data
(optional) ROW_LENGTHS: list of int # the lengths of the sequences in the rows
(optional) COL_LENGTHS: list of int # not used if SYMMETRIC LABELS is True
```

Sparse matrices are stored as canonical CSR arrays: explicit zero values are removed and indices are sorted when written or loaded. Dense files preserve every matrix cell, including zero values, and are usually preferable when most entries are non-zero. Sparse files are usually preferable for sequence similarity graphs or other matrices where most entries are zero.

`MATRIX_FILE_VERSION` is the compatibility marker for this schema. It should be incremented when a change would prevent older Domainator versions from reading a file correctly or would change the meaning of an existing dataset or attribute.

In the future, we may add a `DESCRIPTION` field to distinguish different kinds of data, like raw scores, normalized scores, etc. But so far it is up to the user to remember what the data represents. We may also add `ROW_SEQ_LENGTHS` and `COL_SEQ_LENGTHS` variables to allow for calculation of scores using the EFI score formula, which normalizes on length.

## Tabular data matrix

In addition to hdf5 formatted tabular data, Domainator also supports plain text dense matrices. 
Plain text matrix output is generally specified using the `--dense_text [filename]` command line argument, for example by `seq_dist.py`. For matrix input, the data format is automatically detected.

Domainator plain text matrices are tab-separated and require row and column labels. For example:
```
    seq_1   seq_2
seq_1   1.0 0.3
seq_2   0.3 1.0
```

## Color specification

Some programs, for example, `color_genbank.py` accept color specifications as tab separated text files.

Color files should NOT have a header. 

The columns are `domain_name` and `color`. Where color is a 6-digit RBG hex code, that is not caps-sensitive, and where the `#` prefix is optional.

```python
CcdB	#ff0000
APH	#00ff00
CAT	#0000ff
Condensation	#ff00ff
2-oxoacid_dh	#ffffff
```

## domain_search databases: shards and `.didx` indexes

A *domain_search database* is a GenBank or FASTA file, optionally BGZF-compressed,
optionally split into shards. `domainator_format_db.py` produces these; `domain_search.py`
consumes them. They are created and read transparently — you normally just pass the
logical database name.

### Shard naming

A shard inserts `.N` after the base name, before the format and compression suffixes:

```
mydb.gb         ->  mydb.0.gb,      mydb.1.gb,      ...
mydb.gb.bgz     ->  mydb.0.gb.bgz,  mydb.1.gb.bgz,  ...
```

Leading zeros in the shard number are allowed but ignored (`mydb.00.gb` == shard 0;
they are sorted numerically). `domain_search.py` expands a logical name (`mydb.gb`) to
its shards automatically. If both an unsharded file and shards exist for the same name,
the tools refuse to guess and raise an error — keep only one set.

### `.didx` offset index

An index caches the per-record offsets that `domain_search` would otherwise recompute by
scanning the whole database each run. It is a sidecar named by appending `.didx` to the
*full* file/shard name (`mydb.0.gb.bgz.didx`). Indexes are only supported for uncompressed
or BGZF databases (plain gzip has no usable random-access offset space).

The format is a little-endian binary file: a 56-byte header followed by `record_count`
pairs of two `uint64`.

| offset | field | type | notes |
|--------|-------|------|-------|
| 0  | magic            | `4s`  | `DIDX` |
| 4  | version_major    | `<H`  | current `_INDEX_FILE_VERSION` major (see `db_index.py`) |
| 6  | version_minor    | `<H`  | |
| 8  | filetype_flag    | `<B`  | 0 = genbank, 1 = fasta |
| 9  | compression_flag | `<B`  | 0 = none (byte offsets), 1 = bgzf (virtual offsets `block_start<<16 \| within`) |
| 10 | reserved         | `6x`  | zero pad (body 8-byte aligned) |
| 16 | source_size      | `<Q`  | source file size — staleness fingerprint |
| 24 | source_mtime_ns  | `<Q`  | source mtime (ns) — staleness fingerprint |
| 32 | record_count     | `<Q`  | number of `(offset, cds_count)` pairs in the body |
| 40 | total_cds_count  | `<Q`  | sum of the body's `cds_count` column — lets `domain_search -Z 0` read the target count from the header without scanning the body |
| 48 | flags            | `<Q`  | reserved bitfield |
| 56+| body             | `<QQ` × N | `(offset, cds_count)` pairs, in file order |

A reader **never trusts a stale or unreadable index**: if the magic, version, file kind,
size, or the source size+mtime fingerprint do not match, it warns once and recomputes the
offsets by scanning (an out-of-date offset would otherwise seek into the middle of a
record). Bump the major version in `db_index.py` (`_INDEX_FILE_VERSION`) on any breaking
layout change; a newer minor with the same major stays readable.

### How `domain_search` uses these files

On every run `domain_search` turns each input database into **partitions** —
`(file, offset, n_records)` tuples that worker processes seek to and parse
independently. The offsets come from `utils.get_offsets`/`i_get_offsets`, which read a
valid, fresh `.didx` when one exists and otherwise scan the file (Rust scanner, or
Biopython fallback). Two facts about this drive the layout choices below:

- **A single file already parallelizes across all cores.** Partitions are cut at
  `~--batch_size` CDSs (default 10000), so one large file becomes many partitions spread
  over a pool of `--cpu - 1` workers (the parent does the result merge + GenBank writing).
  Each worker `seek()`s straight to its records — a byte offset for uncompressed files, a
  BGZF virtual offset for compressed ones. **Sharding is not what gives you multi-core
  throughput on one machine; the seek-based partitions already do.**
- **The index removes a serial, parent-side cost.** Without `.didx`, the parent must scan
  each file to produce partitions, and in the default lazy mode that scan runs in the
  single parent thread, which can starve idle workers on a large database. With an index,
  partition production is just a small sidecar read. With `-Z 0` (use the true target
  count for E-values), the count is read straight from the `.didx` headers
  (`total_cds_count`), so an indexed database avoids the otherwise-eager up-front count
  pass entirely.

### Choosing a database layout

**Index (`--index`): almost always, if you will search the database more than once.**
You pay the scan once at build time instead of on every search. Especially worthwhile for
BGZF databases — indexing erases BGZF's only real penalty (its offset *scan* is several
times slower than uncompressed because it must decompress; decompression *during* search
is already at parity). Skip it only for a one-shot search of a database you will discard.

**Compress (`--compress`, BGZF): for anything stored long-term.** ~3–4× smaller on disk
and — unlike plain gzip — still randomly seekable, so parallel search is unaffected.
Combine with `--index` for small files *and* fast startup. Don't bother for tiny or
short-lived databases. (Plain gzip is readable as a convenience but is **not** seekable or
indexable, so it forces slow single-stream reads — re-compress as BGZF for any real use.)

**Shard (`--shards N`): for distribution and operations, not single-machine speed.**
- *One machine:* you generally don't need shards — an indexed single file already
  saturates your cores.
- *Cluster / multi-node:* shard so each node or job owns independent file(s); point each
  `domain_search` at its own shards. Choose the shard count to match your number of
  parallel jobs/nodes, **not** your CPU count.
- *Operational wins:* parallel download/transfer of a large database; replacing or
  re-indexing one slice without touching the whole thing (each shard carries its own
  `.didx`, and a stale shard only invalidates its own index); keeping individual files to a
  manageable size.

Keep only one of {unsharded file, shards} in a directory — `domain_search` resolves a
logical name (`mydb.gb`) to its shards (`mydb.0.gb`, `mydb.1.gb`, …) and refuses to guess
if both are present.

Quick recipes:
- Large reference DB queried repeatedly on a workstation: `--compress --index` (no shards).
- Same DB on an 8-node cluster: `--shards 8 --compress --index`, one shard per node.
- Throwaway one-off search: don't format it at all.