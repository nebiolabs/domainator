# Tree fixtures

These back `test_build_tree.py` and `test_viewer_tree_browser.py`.

## `FeSOD_20.afa`

The 20 sequences of `FeSOD_20.fasta`, aligned to the `Sod_Fe_N` profile of
`FeSOD_pfam.hmm` with pyhmmer's hmmalign, keeping the whole sequences. Their ids are the row
labels of `FeSOD_score_dist.tsv`, so the alignment goes with the tree built from that matrix.

```python
import pyhmmer
from pyhmmer.easel import Alphabet, SequenceFile
from pyhmmer.plan7 import HMMFile

with HMMFile("FeSOD_pfam.hmm") as f:
    hmm = next(h for h in f if h.name == "Sod_Fe_N")
with SequenceFile("FeSOD_20.fasta", digital=True, alphabet=Alphabet.amino()) as f:
    seqs = list(f)
msa = pyhmmer.hmmer.hmmalign(hmm, seqs, trim=False, all_consensus_cols=True)
with open("FeSOD_20.afa", "w") as out:
    for name, aligned in zip(msa.names, msa.alignment):
        aligned = aligned.replace(".", "-").upper()
        out.write(f">{name}\n")
        for i in range(0, len(aligned), 60):
            out.write(aligned[i:i + 60] + "\n")
```

## `tree_support.newick` and `tree_support_metadata.tsv`

A hand-written seven-leaf tree with support values (`90`, `75`), a non-numeric internal label
(`clade_x`) and a multifurcation, and a metadata table for its leaves: a categorical `group`
(with one missing value) and a numeric `size`.
