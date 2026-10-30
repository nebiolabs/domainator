    function treeTopology(view) {
        return Object.fromEntries(['parent', 'length', 'label', 'support', 'leaf_nodes']
            .map(key => [key, view[key].slice()]));
    }

    function validateRootingTopology(view, leaves) {
        const count = view?.parent?.length;
        if (!count || !['parent', 'length', 'label', 'support', 'leaf_nodes'].every(key => Array.isArray(view[key])) ||
                ['length', 'label', 'support'].some(key => view[key].length !== count) ||
                view.parent[0] !== -1 || view.length[0] !== 0 || view.leaf_nodes.length !== leaves) {
            throw new Error('Malformed rooting topology.');
        }
        const children = Array.from({length: count}, () => []);
        const distance = new Float64Array(count);
        const ancestors = [0];
        for (let node = 0; node < count; node++) {
            if (!Number.isFinite(view.length[node]) || view.length[node] < 0) {
                throw new Error('Rooting requires finite, nonnegative branch lengths.');
            }
            if ((view.label[node] !== null && typeof view.label[node] !== 'string') ||
                    (view.support[node] !== null && !Number.isFinite(view.support[node]))) {
                throw new Error('Malformed tree annotations.');
            }
            if (node === 0) { continue; }
            const parent = view.parent[node];
            if (!Number.isInteger(parent) || parent < 0 || parent >= node) {
                throw new Error('Malformed tree parents.');
            }
            while (ancestors.length && ancestors[ancestors.length - 1] !== parent) { ancestors.pop(); }
            if (!ancestors.length) { throw new Error('Tree nodes are not in preorder.'); }
            children[parent].push(node);
            ancestors.push(node);
            distance[node] = distance[parent] + view.length[node];
            if (!Number.isFinite(distance[node])) { throw new Error('Tree distances overflow while rooting.'); }
        }
        const terminals = new Set(children.flatMap((items, node) => items.length ? [] : [node]));
        if (terminals.size !== leaves || new Set(view.leaf_nodes).size !== leaves ||
                view.leaf_nodes.some(node => !terminals.has(node))) {
            throw new Error('Rooting topology does not match the bundle leaves.');
        }
        return children;
    }

    function treeRootGraph(original) {
        const children = validateRootingTopology(original, original.leaf_nodes.length);
        const graph = Array.from({length: original.parent.length}, () => []);
        for (let node = 1; node < graph.length; node++) {
            const edge = {a: original.parent[node], b: node, length: original.length[node], support: original.support[node]};
            graph[edge.a].push({node, edge});
            graph[node].push({node: edge.a, edge});
        }
        // A unary root's stem contains no leaves on its outer side. Leaving it in the
        // unrooted graph would turn the old root into a spurious leaf.
        let start = 0;
        const omitted = [];
        while (children[start].length === 1) {
            omitted.push(start);
            const next = children[start][0];
            graph[next] = graph[next].filter(item => item.node !== start);
            graph[start] = [];
            start = next;
        }
        return {graph, start, omitted};
    }

    function treeMidpointLocation(original) {
        const {graph} = treeRootGraph(original);
        if (original.leaf_nodes.length < 2) { throw new Error('Midpoint rooting needs at least two leaves.'); }
        const leaves = new Set(original.leaf_nodes);
        const walk = start => {
            const parent = new Int32Array(graph.length).fill(-1);
            const distance = new Float64Array(graph.length);
            const stack = [start];
            let farthest = start;
            while (stack.length) {
                const node = stack.pop();
                if (leaves.has(node) && (distance[node] > distance[farthest] ||
                        (distance[node] === distance[farthest] && node < farthest))) {
                    farthest = node;
                }
                for (const item of graph[node]) {
                    if (item.node === parent[node]) { continue; }
                    parent[item.node] = node;
                    distance[item.node] = distance[node] + item.edge.length;
                    if (!Number.isFinite(distance[item.node])) { throw new Error('Tree distances overflow while rooting.'); }
                    stack.push(item.node);
                }
            }
            return {farthest, parent, distance};
        };
        const first = walk(original.leaf_nodes[0]).farthest;
        const {farthest, parent, distance} = walk(first);
        if (!(distance[farthest] > 0)) { throw new Error('Midpoint rooting needs a positive tree diameter.'); }
        const target = distance[farthest] / 2;
        let node = farthest;
        while (parent[node] >= 0 && distance[parent[node]] > target) { node = parent[node]; }
        if (distance[node] === target) { return {node}; }
        const ancestor = parent[node];
        if (distance[ancestor] === target) { return {node: ancestor}; }
        return {a: ancestor, b: node, offset: target - distance[ancestor]};
    }

    function treeRootedTopology(original, location) {
        if (location === null || location?.node === 0) {
            validateRootingTopology(original, original.leaf_nodes.length);
            return {view: treeTopology(original), origins: original.parent.map((_, node) => node), omitted: []};
        }
        const {graph, omitted} = treeRootGraph(original);
        let root;
        if (Number.isInteger(location?.node) && graph[location.node]?.length > 1) {
            root = location.node;
        } else {
            const {a, b, offset} = location || {};
            const item = Number.isInteger(a) && graph[a]?.find(item => item.node === b);
            if (!item || !Number.isFinite(offset) || offset < 0 || offset > item.edge.length) {
                throw new Error('Invalid root position.');
            }
            root = graph.length;
            graph.push([]);
            for (const [node, length] of [[a, offset], [b, item.edge.length - offset]]) {
                graph[node] = graph[node].filter(other => other.node !== (node === a ? b : a));
                const edge = {a: root, b: node, length, support: node === item.edge.b ? item.edge.support : null};
                graph[root].push({node, edge});
                graph[node].push({node: root, edge});
            }
        }
        const view = {parent: [], length: [], label: [], support: [], leaf_nodes: new Array(original.leaf_nodes.length)};
        const bundleIndex = new Map(original.leaf_nodes.map((node, index) => [node, index]));
        const origins = [];
        const stack = [{node: root, from: -1, parent: -1, edge: null}];
        while (stack.length) {
            const {node, from, parent, edge} = stack.pop();
            const index = view.parent.length;
            origins.push(node === original.parent.length ? null : node);
            view.parent.push(parent);
            view.length.push(edge ? edge.length : 0);
            view.label.push(original.label[node] ?? null);
            view.support.push(edge ? edge.support : null);
            if (bundleIndex.has(node)) { view.leaf_nodes[bundleIndex.get(node)] = index; }
            for (let i = graph[node].length - 1; i >= 0; i--) {
                const item = graph[node][i];
                if (item.node !== from) { stack.push({node: item.node, from: node, parent: index, edge: item.edge}); }
            }
        }
        validateRootingTopology(view, original.leaf_nodes.length);
        return {view, origins, omitted};
    }

    function treeSelectedClade(tree, selection, includeRoot = false) {
        if (!selection.size || (!includeRoot && selection.size === tree.view.leaf_nodes.length)) { return null; }
        const selected = new Int32Array(tree.count);
        for (const index of selection) {
            if (!Number.isInteger(index) || index < 0 || index >= tree.view.leaf_nodes.length) { return null; }
            selected[tree.view.leaf_nodes[index]] = 1;
        }
        for (let node = tree.count - 1; node >= 0; node--) {
            if (selected[node] === selection.size && tree.leafCount[node] === selection.size) { return node; }
            if (node > 0) { selected[tree.view.parent[node]] += selected[node]; }
        }
        return null;
    }

    // The original directed edge whose leaf-side is this clade, even when its
    // parent is a synthetic root inserted partway along that edge.
    function treeCladeBoundary(tree, node) {
        const origin = tree.rootOrigins[node];
        const parent = tree.rootOrigins[tree.view.parent[node]];
        if (parent !== null) { return origin + ':' + parent; }
        const location = tree.rootLocation;
        return origin + ':' + (origin === location.a ? location.b : location.a);
    }

    function treeSelectionRootLocation(tree, node) {
        const a = tree.rootOrigins[node];
        const parent = tree.rootOrigins[tree.view.parent[node]];
        if (parent !== null) { return {a, b: parent, offset: tree.view.length[node] / 2}; }
        const location = tree.rootLocation;
        const original = tree.rootOriginal;
        const length = original.parent[location.a] === location.b ? original.length[location.a] : original.length[location.b];
        return {...location, offset: a === location.a ? location.offset / 2 : location.offset / 2 + length / 2};
    }

    function treeRootingProblem(tree) {
        if (tree.view.leaf_nodes.length < 2) { return 'Rooting needs at least two leaves.'; }
        if (tree.view.length.some(length => !Number.isFinite(length) || length < 0)) {
            return 'Rooting requires finite, nonnegative branch lengths.';
        }
        return '';
    }

    function treeMidpointChoice(tree) {
        if (tree.rootMidpoint === null) {
            try {
                tree.rootMidpoint = {location: treeMidpointLocation(tree.rootOriginal), reason: ''};
            } catch (error) {
                tree.rootMidpoint = {location: null, reason: error.message};
            }
        }
        return tree.rootMidpoint;
    }

    function updateTreeRootControls() {
        const tree = state.tree;
        if (!tree || state.activeViewName !== 'tree') { return; }
        const problem = treeRootingProblem(tree);
        const midpoint = document.getElementById('tree-root-midpoint');
        const midpointProblem = problem || treeMidpointChoice(tree).reason;
        midpoint.disabled = Boolean(midpointProblem);
        midpoint.title = midpointProblem || 'Root halfway along the longest leaf-to-leaf path, using branch lengths as given';
        const selection = document.getElementById('tree-root-selection');
        const node = treeSelectedClade(tree, state.selectedNodeIndices);
        selection.disabled = Boolean(problem) || node === null;
        selection.title = problem || (node === null
            ? 'Select exactly one complete clade (or one leaf), but not the entire tree'
            : 'Use the selected clade as the outgroup, rooting halfway along its boundary branch');
        document.getElementById('tree-root-original').disabled = tree.rootLocation === null;
        document.getElementById('tree-root-note').textContent = midpointProblem ||
            (node === null ? 'Selection rooting requires one complete clade or leaf, not the entire tree.' : '');
    }

    function installRootedTree(result, location, message) {
        const previous = state.tree;
        const boundaries = new Set(Array.from(previous.collapsed, node => treeCladeBoundary(previous, node)));
        const childOrders = new Map(Array.from(previous.childOrder, ([node, children]) => [
            node === 0 ? 'root' : treeCladeBoundary(previous, node),
            children.map(child => previous.rootOrigins[child]),
        ]));
        const view = {...previous.view, ...result.view};
        const tree = createTreeState({...state.bundle, views: {tree: view}});
        tree.rootOriginal = previous.rootOriginal;
        tree.rootLocation = location;
        tree.rootOrigins = result.origins;
        tree.rootMidpoint = previous.rootMidpoint;
        tree.rowHeight = previous.rowHeight;
        tree.msaScrollX = previous.msaScrollX;
        for (let node = 0; node < tree.count; node++) {
            const key = node === 0 ? 'root' : treeCladeBoundary(tree, node);
            if (node > 0 && tree.children[node].length && boundaries.has(key)) { tree.collapsed.add(node); }
            const order = childOrders.get(key);
            if (order && order.length === tree.children[node].length) {
                const byOrigin = new Map(tree.children[node].map(child => [tree.rootOrigins[child], child]));
                if (order.every(origin => byOrigin.has(origin))) {
                    tree.childOrder.set(node, order.map(origin => byOrigin.get(origin)));
                }
            }
        }
        if (previous.pendingFrame !== null) { cancelAnimationFrame(previous.pendingFrame); }
        state.bundle.views.tree = view;
        state.tree = tree;
        document.getElementById('tree-cut').value = String(TREE_CUT_DEFAULT);
        document.getElementById('tree-show-support').disabled = !tree.hasSupport;
        treeComputeLayout();
        updateTreeCutLabel();
        renderTreeView();
        updateMetadataTable();
        const expanded = previous.collapsed.size - tree.collapsed.size;
        const clearedOrders = previous.childOrder.size - tree.childOrder.size;
        const omittedAnnotations = result.omitted.some(node => previous.rootOriginal.label[node] ||
            previous.rootOriginal.support[node] !== null) ||
            (location !== null && location.node !== 0 && previous.rootOriginal.support[0] !== null);
        setStatus(message + ' Cluster cut reset to 0 (cut line hidden).' +
            (expanded > 0 ? ' Expanded ' + expanded + ' group(s) that are no longer clades.' : '') +
            (clearedOrders > 0 ? ' Cleared ' + clearedOrders + ' manual child order(s) that no longer apply.' : '') +
            (omittedAnnotations ? ' Root-only annotations are retained in the original-root snapshot, not reassigned.' : ''));
    }

    function rerootTree(mode) {
        if (!state.tree) { return; }
        try {
            const tree = state.tree;
            let location = null;
            if (mode !== 'original') {
                const problem = treeRootingProblem(tree);
                if (problem) { throw new Error(problem); }
                if (mode === 'midpoint') {
                    const choice = treeMidpointChoice(tree);
                    if (choice.reason) { throw new Error(choice.reason); }
                    location = choice.location;
                } else {
                    const node = treeSelectedClade(tree, state.selectedNodeIndices);
                    if (node === null) { throw new Error('Select exactly one complete clade or leaf, but not the entire tree.'); }
                    location = treeSelectionRootLocation(tree, node);
                }
            }
            if (location?.node === 0) { location = null; }
            const result = treeRootedTopology(tree.rootOriginal, location);
            installRootedTree(result, location, mode === 'original' ? 'Restored original root.' :
                mode === 'midpoint' ? 'Rooted at midpoint.' : 'Rooted on selection as outgroup.');
        } catch (error) {
            console.error(error);
            setStatus('Could not reroot tree: ' + error.message);
        }
    }

    function restoreTreeRootingState(value) {
        if (!value || !Object.hasOwn(value, 'location')) { throw new Error('Malformed saved rooting state.'); }
        validateRootingTopology(value.original, state.bundle.nodes.length);
        const result = treeRootedTopology(value.original, value.location);
        for (const key of ['parent', 'length', 'label', 'support', 'leaf_nodes']) {
            if (JSON.stringify(result.view[key]) !== JSON.stringify(state.tree.view[key])) {
                throw new Error('Saved rooting state does not match this tree.');
            }
        }
        state.tree.rootOriginal = treeTopology(value.original);
        state.tree.rootLocation = value.location === null ? null : {...value.location};
        state.tree.rootOrigins = result.origins;
        state.tree.rootMidpoint = null;
    }
