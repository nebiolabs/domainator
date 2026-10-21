"""Checks of the SSN viewer's JavaScript that run in Node, without a browser.

The browser tests (test_ssn_viewer_browser.py) drive the whole page; these check that the
assembled page script parses, and test the layout Web Worker on its own. They skip when
Node is not installed.
"""
import json
import re
import shutil
import subprocess

import pytest

from domainator.ssn_viewer_html import _layout_worker_js, ssn_viewer_html

pytestmark = pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed")

# Runs the worker source with a stand-in `self`, posts it one message, and prints the reply.
_RUN_WORKER = r"""
const fs = require("fs");
const vm = require("vm");
const replies = [];
const context = {self: {postMessage: message => replies.push(message)}, Math, Map, Set, Array, Object, Number, Infinity};
vm.createContext(context);
vm.runInContext(fs.readFileSync(process.argv[2], "utf8"), context);
context.self.onmessage({data: JSON.parse(fs.readFileSync(process.argv[3], "utf8"))});
process.stdout.write(JSON.stringify(replies));
"""


def _node(tmp_path, *args):
    return subprocess.run(["node", *args], capture_output=True, text=True, check=True, cwd=tmp_path).stdout


def test_page_script_and_worker_parse(tmp_path):
    page = ssn_viewer_html()
    scripts = re.findall(r"<script>\n(.*?)\n</script>", page, re.S)
    assert len(scripts) == 1
    for name, source in (("page.js", scripts[0]), ("worker.js", _layout_worker_js())):
        (tmp_path / name).write_text(source)
        _node(tmp_path, "--check", name)


def _layout(tmp_path, visible_ids, links):
    (tmp_path / "worker.js").write_text(_layout_worker_js())
    (tmp_path / "run.js").write_text(_RUN_WORKER)
    message = {
        "type": "computeLayout", "requestId": 7, "key": "k", "algorithm": "force", "visibleIds": visible_ids,
        "links": [{"sourceId": s, "targetId": t} for s, t in links],
        "hierarchyNodes": {str(i): {"size": 1, "leaf_start": i} for i in visible_ids},
        "sortBySizeEnabled": True,
    }
    (tmp_path / "message.json").write_text(json.dumps(message))
    (reply,) = json.loads(_node(tmp_path, "run.js", "worker.js", "message.json"))
    assert (reply["requestId"], reply["key"]) == (7, "k")
    return {item["componentId"]: item for item in reply["layout"]}


def test_worker_lays_out_every_node_without_overlaps(tmp_path):
    # A chain of four, a pair, and a singleton: three connected components.
    ids = list(range(7))
    layout = _layout(tmp_path, ids, [(0, 1), (1, 2), (2, 3), (4, 5)])
    assert sorted(layout) == ids
    for item in layout.values():
        assert all(isinstance(item[key], (int, float)) and abs(item[key]) < 1e6 for key in ("x", "y", "radius"))
    for a in ids:
        for b in ids[a + 1:]:
            distance = ((layout[a]["x"] - layout[b]["x"]) ** 2 + (layout[a]["y"] - layout[b]["y"]) ** 2) ** 0.5
            assert distance >= layout[a]["radius"] + layout[b]["radius"] - 1e-6, (a, b)


def test_worker_layout_is_deterministic(tmp_path):
    ids = list(range(6))
    links = [(0, 1), (1, 2), (0, 2), (3, 4)]
    assert _layout(tmp_path, ids, links) == _layout(tmp_path, ids, links)
