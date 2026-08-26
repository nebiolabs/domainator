# Structure test fixtures

Small structures for the `structure_*` tools. Each file was downloaded from RCSB and
trimmed to keep the repository small:

* first model only (matters for NMR entries),
* altloc blank or `A` only,
* backbone + CB atoms only (`N`, `CA`, `C`, `O`, `CB`),
* no `HETATM`, waters, ligands or hydrogens.

foldseek derives its 3Di states from virtual centers built on `N`/`CA`/`CB`, so
backbone + CB is the smallest atom set that exercises the normal code path. (C-alpha only
also works and produces near-identical 3Di, if these ever need to shrink further.)

Regenerate with:

```bash
curl -O https://files.rcsb.org/download/<ID>.pdb
python3 scripts/trim_pdb_fixtures.py <ID>.pdb <ID>_bb.pdb N,CA,C,O,CB
gzip -9 <ID>_bb.pdb
```

## Contents

| file | PDB | chains | residues | role |
| ---- | --- | ------ | -------- | ---- |
| `refs/UBQref.pdb.gz` | 1UBQ | A | 76 | the reference structure searched with |
| `inputs/1UBQ.pdb.gz` | 1UBQ | A | 76 | identical to the reference: 100% identity self-hit |
| `inputs/1NDD.pdb.gz` | 1NDD | A | 74 | NEDD8, a ubiquitin homolog: a real remote hit (~57% identity) |
| `inputs/1CRN.pdb.gz` | 1CRN | A | 46 | crambin, unrelated fold: no hit |
| `inputs/1IGD.pdb.gz` | 1IGD | A | 61 | protein G B3 domain, unrelated fold: no hit |
| `inputs/1ZNI.pdb.gz` | 1ZNI | A-D | 21/30/21/30 | insulin, multi-chain: one record per chain, no hit |
| `inputs/some.long.name.pdb.gz` | 1IGD | A | 61 | copy of 1IGD under a dotted filename, to pin entry naming and GenBank LOCUS handling |

`inputs/` is deliberately searchable as a directory: `1UBQ` and `1NDD` hit `UBQref` and
the rest do not, so a single search exercises both the hit and no-hit paths.
