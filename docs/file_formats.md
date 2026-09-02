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


Annotations produced by a structural aligner (`structure_domainate.py`, `structure_search.py`) carry the same qualifiers, with `program` naming the aligner (for example `foldseek`). They may additionally carry `tmscore`, `lddt`, `rmsd` and `prob` when those were requested with `--metrics`. These four qualifiers are optional and are simply absent when the aligner did not compute them, so files written by the sequence-based tools are unchanged. They are scores for the alignment as a whole, on their own scales rather than on the `score` bitscore scale, so they should not be compared against `score` or against each other.

```javascript
     Domainator      1..76
                     /program="foldseek"
                     /database="refs"
                     /description="."
                     /evalue="3.1e-18"
                     /score="635.0"
                     /name="UBQref_A"
                     /identity="100.0"
                     /cds_id="0_1_76"
                     /rstart="1"
                     /rend="76"
                     /rlen="76"
                     /tmscore="1.013"
                     /lddt="1.000"
                     /rmsd="0.015"
                     /prob="1.000"
```

Note that for structure-derived records, all coordinates are indices into the *observed* residue sequence that the aligner extracted from the structure file, not original PDB residue numbers. A chain with a disordered loop has a sequence shorter than its residue-number span, so positions cannot be mapped back to the structure's own numbering.

`misc_feature` features carrying a `threedi` qualifier are written by `structure_to_genbank.py --store_3di`. The qualifier holds the chain's full-length 3Di structural-alphabet string, one character per residue. It is a single opaque string rather than a per-residue annotation, so it does not survive tools that slice records (`extract_domains.py`, `select_by_coord.py`).

`Domain_Search` features are added by `domain_search.py` and `structure_search.py`. `Domain_Search` annotations are related to `Domainator` annotations in that they have exactly the same qualifiers. What makes them distinct is that they are cleared from sequences used as input to `domain_search.py`, and there are exactly one per contig on sequences returned by `domain_search.py`. Other programs, such as `extract_domains.py`, `select_by_cds.py`, and others also handle `Domain_Search` annotations distinctly from `Domainator` annotations, typically by using the `--search_hits` command line argument.

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

## Protein structures

The `structure_*` tools read protein structure files: `.pdb`, `.pdb1`, `.ent`, `.cif`,
`.mmcif` and `.bcif`, optionally gzipped. Domainator does not parse these itself — the
structural aligner does, which is also what splits multi-chain structures into one record
per chain and names them `<file basename>_<chain>`.

These extensions are recognized by `domainator.utils.get_file_type`, which reports them as
type `structure`. That type is deliberately not readable by any sequence reader: every
consumer of `get_file_type` whitelists the types it accepts, so a structure file handed to
`domainate.py`, `domain_search.py` or any other sequence tool is rejected with a message
pointing at [structure_to_genbank.py and friends](structure_search.md) rather than a generic
"extension not recognized".

Aligner databases are *not* identified by extension. A foldseek database is a set of files
sharing a prefix, detected by the presence of `<prefix>.dbtype` and `<prefix>.index`; a
reseek database is a single `.bcb`/`.bca` file, detected by its magic bytes. Extension-based
detection is avoided for both so that a mis-named file is not mistaken for a database.

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

Sparse matrices are stored as canonical CSR arrays: explicit zero values are removed and indices are sorted when written or loaded. `SPARSE_CSR_INDICES` and `SPARSE_CSR_INDPTR` are written at the narrowest width that can address the matrix (int32 unless the number of rows, columns, or non-zeros exceeds the int32 maximum), so files written by different tools may use different index widths; readers must accept either. Dense files preserve every matrix cell, including zero values, and are usually preferable when most entries are non-zero. Sparse files are usually preferable for sequence similarity graphs or other matrices where most entries are zero.

Threshold tables derived from a matrix (`matrix_report.py`, `build_ssn_viewer.py`) use
one row per distinct maximum-spanning-tree edge weight, and every row describes the graph
that `build_ssn --lb <threshold>` keeps: edges scoring **strictly above** the threshold.
An extra `threshold = 0` row (the `--lb 0` default) closes the table out with the complete
graph when every score is positive. Ties at a threshold are excluded from that
threshold's own row, so counts are a function of the threshold alone rather than of the
order in which equal-weight edges happened to be sorted.

`MATRIX_FILE_VERSION` is the compatibility marker for this schema. It should be incremented when a change would prevent older Domainator versions from reading a file correctly or would change the meaning of an existing dataset or attribute.

In the future, we may add a `DESCRIPTION` field to distinguish different kinds of data, like raw scores, normalized scores, etc. But so far it is up to the user to remember what the data represents. We may also add `ROW_SEQ_LENGTHS` and `COL_SEQ_LENGTHS` variables to allow for calculation of scores using the EFI score formula, which normalizes on length.

## SSN viewer bundles (`.ssnv`)

`build_ssn_viewer.py` turns a symmetric similarity matrix into an `.ssnv` bundle: a
**gzip-compressed JSON document** read by the standalone HTML viewer
(`--html`/`--embed_data`) and by `ssn_navigator.py`. Readers accept plain (uncompressed)
JSON as well, detected by the gzip magic number. The bundle stores an MST-derived merge
hierarchy rather than the matrix itself, so its size scales with the node count, not with
O(n²) edges.

Top-level keys:

```python
format: "domainator_ssn_viewer_bundle"  # a mismatch here is fatal for every reader
version: int                            # see the compatibility rules below
name: str                               # display name
domainator_version: str                 # the build that wrote the file
graph: dict                             # nodes, mst_edges, hierarchy, slider_stops, merge series
metadata: dict                          # positional node metadata table (see below)
defaults: dict                          # color_by, label_by, categorical_columns
(v4, optional) app_state: dict          # viewer UI state -- see below
```

`graph.slider_stops` closes with a floor stop strictly below the weakest MST edge. Every
other stop excludes its own tie group under the strictly-above rule, so without that last
one the fully merged network (the true connected components) could not be reached: the
lowest stop would still split the weakest merge apart. The floor sits **1% of the MST
weight range** below the weakest edge rather than at `0`, so that a network scoring
350–650 does not spend most of its slider track on empty space below the data, and so
that negative scores are cleared too. It reports the `threshold_index` of the
complete-graph row, since that is the cut it stands for.

`graph.nodes` is the list of node ids, and its **position is the node index** used
throughout the bundle: `metadata.rows[i]` describes `graph.nodes[i]`, and
`graph.hierarchy.leaf_order` holds those same indices. `metadata` is a positional table
rather than a per-row mapping:

```python
metadata: {
  "columns": [{"name": str, "type": "str" | "int" | "float", (optional) "origin": str}, ...],
  "rows":    [[value, value, ...], ...]   # one row per node, in graph.nodes order
}
```

`origin: "viewer"` marks a column created in the HTML viewer rather than merged from a
`--metadata` TSV; it is what lets the viewer offer to delete the column again. Readers
that do not know the key can ignore it.

### Versions and compatibility

The version constants live in `ssn_bundle.py`
(`SSN_VIEWER_BUNDLE_VERSION`, `SUPPORTED_SSN_VIEWER_BUNDLE_VERSIONS`).

| Version | Change |
| --- | --- |
| 3 | Per-event `merge_size_counts`/`largest_merge`/`merge_count` and `graph.merge_moving_sum`. |
| 4 | Optional top-level `app_state`. Purely additive: `build_ssn_viewer.py` never writes it, and a reader that ignores the section can treat a v4 file exactly like a v3 file. |

Bump `SSN_VIEWER_BUNDLE_VERSION` for any change to the schema, and add the old version to
`SUPPORTED_SSN_VIEWER_BUNDLE_VERSIONS` when the change is additive so previously written
bundles keep loading. `ssn_bundle.load_bundle` refuses a version outside that tuple; the
HTML viewer is deliberately more permissive and loads an unknown version with a warning,
because refusing a file it can very likely still render is the worse failure.

### `app_state` (saved sessions)

The viewer's "Save session…" button writes the bundle back out with an `app_state`
section carrying the UI state. **The rest of the file is an ordinary bundle**: metadata
edited in the viewer is written into `metadata.rows`/`metadata.columns`, so
`ssn_navigator.py` and any other reader see the annotations without understanding
`app_state` at all.

```python
app_state: {
  "state_version": int, "saved_by": str, "saved_at": str,   # ISO 8601
  "view":      {...},   # layout, color_by/label_by, toggles, threshold_value, view_transform
  "table":     {...},   # sort, filter, null_order, rows_per_page, column_widths
  "colors":    {...},   # custom_palettes, categorical_columns
  "selection": {"node_ids": [...], "presets": {"<slot 0-9>": {"node_ids": [...], "saved_at": str}}}
}
```

The viewer can also write an **extraction**: a bundle over a chosen subset of nodes,
rebuilt from the induced MST edges rather than re-derived from the source matrix.
`ssn_viewer_html.py` ports `ssn_hierarchy.py`'s hierarchy and merge-series maths into
JavaScript to do this, so the two must be changed together; the browser tests assert the
port reproduces `build_ssn_viewer.py`'s output field-for-field when the whole network is
selected.

An extraction is only written when every original network component contributes a single
MST-connected piece to the selection. The MST kept one path between any two nodes and
discarded the rest, so a selection that omits the nodes along that path leaves the bundle
with no evidence of how the surviving pieces relate — writing them out as separate
components would assert an absence of similarity the file cannot support. Nodes in
different components of the original network are exempt, since they were already
unrelated. To extract an arbitrary subset, subset the source matrix instead
(`build_ssn_viewer.py --subset`), which measures those relationships rather than
inferring them. One field cannot be rebuilt this way and is written empty:
`graph.edges_by_threshold` counts edges of the *full* graph, which a bundle never carried.

Two rules make this section survive UI churn, and changes to it should preserve both:

- **Every field is optional and independently skippable.** The viewer defines one registry
  entry (`VIEW_STATE_FIELDS` in `ssn_viewer_html.py`) per persisted control and looks each
  saved key up in it. A key the current build does not recognize, or one whose value it
  cannot apply, is skipped and reported in the load status rather than aborting the load; a
  key the file omits keeps its default. So a session written by a different version of the
  viewer still opens.
- **Selections are stored as node ids, never indices.** An index is a position in
  `graph.nodes` and would silently mean a different sequence if the bundle were rebuilt or
  subsetted. Ids that are not in the bundle are counted and reported on load. For the same
  reason the threshold is stored as `threshold_value`, not as a slider position, which is
  derived from `graph.slider_stops` and shifts whenever `--max_merge_events` changes.


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