/*
 * The layout Web Worker's entry point. The worker's source is layout_core.js
 * followed by this file (see viewer_html._layout_worker_js).
 */

self.onmessage = function(event) {
    const msg = event.data;
    if (msg.type !== 'computeLayout') { return; }
    const {requestId, key, algorithm, visibleIds, links, hierarchyNodes, sortBySizeEnabled} = msg;
    const {adjacency, components} = clusterGraphComponents(visibleIds, links, hierarchyNodes, sortBySizeEnabled);
    const componentLayouts = components.map(ids =>
        simulateComponentLayout(ids, adjacency, forestForceOptions(ids.length), hierarchyNodes));
    const layout = packLayouts(componentLayouts, FOREST_PACK_OPTIONS);
    self.postMessage({requestId, key, layout});
};
