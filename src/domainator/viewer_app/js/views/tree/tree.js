    // ---- Tree view: the tree and its layout ----
    //
    // A bundle's views.tree holds the tree as arrays over its nodes in preorder (parent, length,
    // label, support), so a node's subtree is the run of nodes from it up to its preorder end,
    // and leaf_nodes[i] is the tree node of the bundle's node i. The layout is rectangular: one
    // row per visible leaf (a collapsed clade takes one row), each internal node halfway between
    // its first and last visible child, and x the distance from the root (a phylogram) or, in a
    // cladogram, the number of steps above the deepest leaf.

    const TREE_ROW_HEIGHT_DEFAULT = 16;
    const TREE_ROW_HEIGHT_MAX = 40;
    const TREE_ROW_HEIGHT_MIN = 0.01;
    // The slider's positions; its last position means no cut.
    const TREE_CUT_STEPS = 1000;

    function createTreeState(bundle) {
        const view = bundle.views.tree;
        const count = view.parent.length;
        const children = Array.from({length: count}, () => []);
        for (let node = 1; node < count; node++) {
            children[view.parent[node]].push(node);
        }
        const treeToBundle = new Int32Array(count).fill(-1);
        view.leaf_nodes.forEach((treeNode, bundleIndex) => { treeToBundle[treeNode] = bundleIndex; });
        const subtreeEnd = new Int32Array(count);
        const leafCount = new Int32Array(count);
        const rootDistance = new Float64Array(count);
        const height = new Int32Array(count);
        for (let node = 1; node < count; node++) {
            rootDistance[node] = rootDistance[view.parent[node]] + Math.max(0, Number(view.length[node]) || 0);
        }
        for (let node = count - 1; node >= 0; node--) {
            if (children[node].length === 0) {
                leafCount[node] = 1;
                subtreeEnd[node] = node + 1;
            } else {
                const last = children[node][children[node].length - 1];
                subtreeEnd[node] = subtreeEnd[last];
                children[node].forEach(child => {
                    leafCount[node] += leafCount[child];
                    height[node] = Math.max(height[node], height[child] + 1);
                });
            }
        }
        return {
            view, count, children, treeToBundle, subtreeEnd, leafCount, rootDistance, height,
            hasSupport: view.support.some(value => value !== null && value !== undefined) ||
                view.label.some(label => label !== null && label !== undefined && label !== ''),
            msa: view.msa || null,
            collapsed: new Set(),
            // Pixels per row, and how far the rows are scrolled; null until the first paint fits
            // the tree to the canvas.
            rowHeight: null,
            scrollY: 0,
            msaScrollX: 0,
            layout: null,
            msaConsensus: null,
            cladeConsensus: new Map(),
            pendingFrame: null,
            // The selection box being dragged: {canvas, x0, y0, x1, y1}.
            box: null,
        };
    }

    function treeIsLeaf(node) {
        return state.tree.children[node].length === 0;
    }

    // The bundle's nodes (leaves) in the subtree of `node`.
    function treeCladeLeaves(node) {
        const tree = state.tree;
        const leaves = [];
        for (let other = node; other < tree.subtreeEnd[node]; other++) {
            if (tree.treeToBundle[other] >= 0) {
                leaves.push(tree.treeToBundle[other]);
            }
        }
        return leaves;
    }

    function treeLayoutMode() {
        return document.getElementById('tree-layout').value;
    }

    function treeLadderized() {
        return document.getElementById('tree-ladderize').checked;
    }

    function treeComputeLayout() {
        const tree = state.tree;
        const count = tree.count;
        const displayChildren = tree.children.map(children => children.slice());
        if (treeLadderized()) {
            displayChildren.forEach(children => children.sort((left, right) => tree.leafCount[left] - tree.leafCount[right]));
        }
        const x = new Float64Array(count);
        if (treeLayoutMode() === 'cladogram') {
            const top = tree.height[0];
            for (let node = 0; node < count; node++) {
                x[node] = top - tree.height[node];
            }
        } else {
            x.set(tree.rootDistance);
        }
        let maxX = 0;
        for (let node = 0; node < count; node++) {
            maxX = Math.max(maxX, x[node]);
        }
        // The farthest x in each clade, where a collapsed clade's wedge ends.
        const cladeMaxX = Float64Array.from(x);
        for (let node = count - 1; node > 0; node--) {
            const parent = tree.view.parent[node];
            cladeMaxX[parent] = Math.max(cladeMaxX[parent], cladeMaxX[node]);
        }

        const rowNode = [];
        const nodeRow = new Float64Array(count).fill(NaN);
        const visibleOrder = [];
        const stack = [0];
        while (stack.length > 0) {
            const node = stack.pop();
            visibleOrder.push(node);
            if (displayChildren[node].length === 0 || tree.collapsed.has(node)) {
                nodeRow[node] = rowNode.length;
                rowNode.push(node);
                continue;
            }
            for (let i = displayChildren[node].length - 1; i >= 0; i--) {
                stack.push(displayChildren[node][i]);
            }
        }
        for (let i = visibleOrder.length - 1; i >= 0; i--) {
            const node = visibleOrder[i];
            if (Number.isNaN(nodeRow[node])) {
                const children = displayChildren[node];
                nodeRow[node] = (nodeRow[children[0]] + nodeRow[children[children.length - 1]]) / 2;
            }
        }
        tree.layout = {x, maxX: maxX || 1, cladeMaxX, rowNode, nodeRow, visibleOrder, displayChildren};
        tree.cladeConsensus = new Map();
    }

    // ---- Cutting the tree into clusters ----

    // The distance from the root the cut is at, in the layout's x units; Infinity when the
    // slider is at its end (no cut).
    function treeCutX() {
        const position = Number(document.getElementById('tree-cut').value);
        if (!state.tree || !state.tree.layout || position >= TREE_CUT_STEPS) {
            return Infinity;
        }
        return position / TREE_CUT_STEPS * state.tree.layout.maxX;
    }

    function updateTreeCutLabel() {
        const cutX = treeCutX();
        document.getElementById('tree-cut-label').textContent =
            Number.isFinite(cutX) ? formatTreeLength(cutX) : '(no cut)';
    }

    function formatTreeLength(value) {
        return Number(value.toPrecision(3)).toString();
    }

    // Every branch that crosses the cut starts a cluster: the clade below it. A leaf left of the
    // cut is a cluster of its own. Clusters are numbered from 1, largest first, as the SSN
    // view numbers its clusters.
    function treeClusterNumbers() {
        const tree = state.tree;
        const cutX = treeCutX();
        const x = tree.layout.x;
        const clusterRoot = new Int32Array(tree.count).fill(-1);
        for (let node = 0; node < tree.count; node++) {
            const parent = tree.view.parent[node];
            if (parent >= 0 && clusterRoot[parent] >= 0) {
                clusterRoot[node] = clusterRoot[parent];
            } else if (x[node] >= cutX || tree.children[node].length === 0) {
                clusterRoot[node] = node;
            }
        }
        const roots = [];
        for (let node = 0; node < tree.count; node++) {
            if (clusterRoot[node] === node) {
                roots.push(node);
            }
        }
        // Preorder index breaks ties, so equal clusters are numbered top to bottom.
        roots.sort((left, right) => tree.leafCount[right] - tree.leafCount[left] || left - right);
        const numberByRoot = new Map(roots.map((root, index) => [root, index + 1]));
        return {
            clusterCount: roots.length,
            numberFor: bundleIndex => numberByRoot.get(clusterRoot[tree.view.leaf_nodes[bundleIndex]]) ?? null,
        };
    }

    // ---- Newick ----

    function newickName(name) {
        return String(name).replace(/[ ;:,()']/g, '_');
    }

    function formatNewickNumber(value) {
        return Number(Number(value).toPrecision(6)).toString();
    }

    // The tree in Newick format, pruned to the bundle nodes in `keep` (all of them when null):
    // other leaves are removed, then internal nodes left with no leaves, and an internal node
    // left with one child is merged into it, adding their branch lengths, as build_tree --subset
    // does.
    function treeNewick(keep) {
        const tree = state.tree;
        const view = tree.view;
        const kept = new Uint8Array(tree.count);
        for (let node = tree.count - 1; node >= 0; node--) {
            if (tree.children[node].length === 0) {
                kept[node] = keep === null || keep.has(tree.treeToBundle[node]) ? 1 : 0;
            } else {
                kept[node] = tree.children[node].some(child => kept[child]) ? 1 : 0;
            }
        }
        if (!kept[0]) {
            return null;
        }
        const label = node => {
            if (tree.children[node].length === 0) {
                return newickName(nodeId(tree.treeToBundle[node]));
            }
            if (view.label[node] !== null && view.label[node] !== undefined && view.label[node] !== '') {
                return newickName(view.label[node]);
            }
            return view.support[node] === null || view.support[node] === undefined ? '' : formatNewickNumber(view.support[node]);
        };
        // Each item is a node to write, with the length carried down from merged ancestors, or
        // a string to write as is.
        const out = [];
        const stack = [{node: 0, carried: 0, root: true}];
        while (stack.length > 0) {
            const item = stack.pop();
            if (typeof item === 'string') {
                out.push(item);
                continue;
            }
            const {node, root} = item;
            const length = item.carried + (node === 0 ? 0 : Number(view.length[node]) || 0);
            const children = tree.children[node].filter(child => kept[child]);
            if (children.length === 1) {
                stack.push({node: children[0], carried: length, root});
                continue;
            }
            const branch = root ? '' : ':' + formatNewickNumber(length);
            if (children.length === 0) {
                out.push(root ? '(' + label(node) + ':0)' : label(node) + branch);
                continue;
            }
            out.push('(');
            stack.push(')' + label(node) + branch);
            for (let i = children.length - 1; i >= 0; i--) {
                stack.push({node: children[i], carried: 0, root: false});
                if (i > 0) {
                    stack.push(',');
                }
            }
        }
        return out.join('') + ';';
    }
