    // ---- Tree view: the alignment ----
    //
    // views.tree.msa is {alphabet, length, sequences}: one aligned sequence per bundle node, or
    // null for a leaf the alignment lacks. The alignment is drawn beside the tree, one row per
    // tree row, so a collapsed clade shows its consensus.

    const MSA_GAPS = new Set(['-', '.']);
    // Amino acids by class, and nucleotides, as in the ancestral-reconstruction viewer.
    const MSA_PROTEIN_COLORS = {
        A: '#96c8f0', I: '#96c8f0', L: '#96c8f0', M: '#96c8f0', F: '#96c8f0', W: '#96c8f0', V: '#96c8f0',
        C: '#f2a0c0', K: '#f0827a', R: '#f0827a', D: '#ea6ee0', E: '#ea6ee0',
        N: '#a0f0a0', Q: '#a0f0a0', S: '#a0f0a0', T: '#a0f0a0',
        G: '#f5b96a', P: '#f5e06a', H: '#7ad8d8', Y: '#7ad8d8',
    };
    const MSA_NUCLEOTIDE_COLORS = {A: '#a0e0a0', C: '#96c8f0', G: '#f5c96a', T: '#f0928a', U: '#f0928a'};
    const MSA_GAP_COLOR = '#f4f2ee';
    const MSA_OTHER_COLOR = '#dedad2';
    const MSA_IDENTICAL_COLOR = '#8fb3cc';
    const MSA_PLAIN_COLOR = '#e3e6ea';

    function msaScheme() {
        const scheme = document.getElementById('msa-scheme').value;
        if (scheme === 'auto') {
            return state.tree.msa.alphabet === 'nucleotide' ? 'nucleotide' : 'protein';
        }
        return scheme;
    }

    // The color of residue `letter` in column `column`, in `scheme`.
    function msaResidueColor(letter, column, scheme) {
        if (MSA_GAPS.has(letter)) {
            return MSA_GAP_COLOR;
        }
        const upper = letter.toUpperCase();
        if (scheme === 'protein') {
            return MSA_PROTEIN_COLORS[upper] || MSA_OTHER_COLOR;
        }
        if (scheme === 'nucleotide') {
            return MSA_NUCLEOTIDE_COLORS[upper] || MSA_OTHER_COLOR;
        }
        if (scheme === 'identity') {
            return upper === msaConsensus().letters[column] ? MSA_IDENTICAL_COLOR : '#ffffff';
        }
        return MSA_PLAIN_COLOR;
    }

    // The most common residue of each column among `sequences` (a gap when a column has none),
    // and the fraction of the column's residues that are it.
    function msaColumnConsensus(sequences, length) {
        const letters = new Array(length).fill('-');
        const identity = new Float64Array(length);
        for (let column = 0; column < length; column++) {
            const counts = new Map();
            let residues = 0;
            sequences.forEach(sequence => {
                const letter = sequence[column].toUpperCase();
                if (!MSA_GAPS.has(letter)) {
                    counts.set(letter, (counts.get(letter) || 0) + 1);
                    residues += 1;
                }
            });
            let best = 0;
            counts.forEach((count, letter) => {
                if (count > best) {
                    best = count;
                    letters[column] = letter;
                }
            });
            identity[column] = residues > 0 ? best / sequences.length : 0;
        }
        return {letters, identity};
    }

    function msaConsensus() {
        const tree = state.tree;
        if (!tree.msaConsensus) {
            const sequences = tree.msa.sequences.filter(sequence => sequence !== null);
            tree.msaConsensus = msaColumnConsensus(sequences, tree.msa.length);
        }
        return tree.msaConsensus;
    }

    // The sequence a row shows: a leaf's own, or a collapsed clade's consensus.
    function msaRowSequence(node) {
        const tree = state.tree;
        if (tree.children[node].length === 0) {
            return tree.msa.sequences[tree.treeToBundle[node]];
        }
        if (!tree.cladeConsensus.has(node)) {
            const sequences = treeCladeLeaves(node).map(index => tree.msa.sequences[index]).filter(sequence => sequence !== null);
            tree.cladeConsensus.set(node, sequences.length > 0 ? msaColumnConsensus(sequences, tree.msa.length).letters.join('') : null);
        }
        return tree.cladeConsensus.get(node);
    }

    // The selected leaves' sequences (every leaf's when none is selected) as FASTA, with or
    // without the alignment's gaps.
    function msaFasta(ungapped) {
        const tree = state.tree;
        const selected = state.selectedNodeIndices.size > 0;
        const lines = [];
        tree.msa.sequences.forEach((sequence, bundleIndex) => {
            if (sequence === null || (selected && !state.selectedNodeIndices.has(bundleIndex))) {
                return;
            }
            const text = ungapped ? sequence.replace(/[-.]/g, '') : sequence;
            lines.push('>' + nodeId(bundleIndex));
            for (let start = 0; start < text.length; start += 60) {
                lines.push(text.slice(start, start + 60));
            }
        });
        return lines.length > 0 ? lines.join('\n') + '\n' : '';
    }
