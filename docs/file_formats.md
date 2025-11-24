[index](README.md)
# Domainator file formats

## Genbank

Genbank files are the sequence format most used by Domainator. Genbank files are convenient because they can store information on both sequence and sequence positional annotations.

### Features
Domainator produces and understands some feature types particular to it.

`Domainator` features are domain annotations added to genbank files by `domainate.py` or `domain_search.py`. The `cds_id` qualifier is particularly critical as it links the `Domainator` feature to a `CDS` feature. Domainator adds a `cds_id` qualifier to any `CDS` feature lacking one. The `cds_id` is expected to be unique within a contig. The general format for an automatically generated `cds_id` is `[start coordinate]_[strand]_[end coordinate]`. Multi-part CDSs, for example CDSs containing introns or crossing the origin are handled in a somewhat more complicated way, as special cases, with the goal being to maintain a one-cds_id to one_cds relation.

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

Domainator uses hdf5 files to store data matrics in both dense and sparse formats. Code for reading, writing, and using data matrices is found in the `data_matrix.py` file (which is not executable, but imported by various other scripts).

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
SPARSE_CSR_INDICES: list of int # row index for CSR sparse data
SPARSE_CSR_INDPTR: list of int # col index for CSR sparse data
(optional) ROW_LENGTHS: list of int # the lengths of the sequences in the rows
(optional) COL_LENGTHS: list of int # not used if SYMMETRIC LABELS is True
```

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