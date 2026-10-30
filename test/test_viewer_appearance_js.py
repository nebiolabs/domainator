"""Pure viewer shape and appearance helpers, exercised without a browser."""

import json
import shutil
import subprocess

import pytest

from domainator.viewer_html import _script

pytestmark = pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed")


def _run(tmp_path, source):
    script = tmp_path / "appearance.js"
    script.write_text(_script("common/shapes.js") + _script("common/appearance.js") + source)
    return json.loads(subprocess.run(["node", str(script)], capture_output=True, text=True, check=True).stdout)


def test_glyph_geometry_and_svg_share_bounded_vertices(tmp_path):
    result = _run(tmp_path, """
        const results = NODE_SHAPES.map(([shape]) => {
            const points = [];
            const ctx = {moveTo: (...p) => points.push(p), lineTo: (...p) => points.push(p),
                closePath() {}, arc: (x, y, r) => points.push([x, y, r])};
            appendNodeGlyph(ctx, shape, 10, 20, 4);
            return {shape, vertices: shapeVertices(shape, 10, 20, 4), points,
                svg: nodeGlyphSVG(shape, 10, 20, 4)};
        });
        process.stdout.write(JSON.stringify(results));
    """)
    assert [row["shape"] for row in result] == [
        "circle", "square", "diamond", "triangle_up", "triangle_down", "hexagon",
    ]
    for row in result:
        if row["shape"] == "circle":
            assert row["vertices"] is None and row["svg"].startswith("<circle")
        else:
            assert row["vertices"] == row["points"] and row["svg"].startswith("<polygon")
            for x, y in row["vertices"]:
                assert (x - 10) ** 2 + (y - 20) ** 2 == pytest.approx(16)


def test_shape_assignment_numeric_categories_nulls_and_overrides(tmp_path):
    result = _run(tmp_path, """
        const values = [10, 2, null, '', 3, 4, 5, 6, 7, '__proto__'];
        const state = {allNodeIndices: values.map((_, i) => i), shapeCaches: new Map(),
            customShapes: {score: {values: {'2': 'hexagon'}, null_shape: 'diamond'}}};
        const sortPaletteKeys = keys => keys.sort((a, b) => Number(a) - Number(b));
        const distinctColumnValues = () => ({values: values.filter(v => v !== null && v !== '').map(v => ({key: String(v)}))});
        const metadataValue = i => values[i];
        const automatic = automaticShapes('score');
        const shapes = values.map((_, i) => shapeForNode({shape_by: 'score'}, i));
        process.stdout.write(JSON.stringify({automatic, shapes,
            fallback: shapeForNode({shape_by: ''}, 0, 'cell')}));
    """)
    assert result["shapes"][1] == "hexagon"
    assert result["shapes"][2:4] == ["diamond", "diamond"]
    assert result["fallback"] == "cell"
    assert result["automatic"]["2"] == "circle"
    assert result["automatic"]["7"] == "hexagon"
    assert result["automatic"]["10"] == "circle"
    assert result["automatic"]["__proto__"] == "square"


def test_style_validation_preserves_empty_tree_lists_and_rejects_bad_data(tmp_path):
    result = _run(tmp_path, """
        const state = {activeViewName: 'tree'};
        const metadataColumn = name => name === 'group';
        const results = [validateAppearanceStyles([], 'node')];
        for (const [view, value] of [
            ['ssn', []],
            ['tree', [{id: 'node-1', color_by: 'gone', shape_by: ''}]],
            ['tree', [{id: 'node-1', color_by: '', shape_by: ''}, {id: 'node-1', color_by: '', shape_by: ''}]],
            ['tree', [{id: 'label-1', color_by: '', shape_by: ''}]]
        ]) {
            state.activeViewName = view;
            try { validateAppearanceStyles(value, 'node'); results.push('accepted'); }
            catch (error) { results.push(error.message); }
        }
        process.stdout.write(JSON.stringify(results));
    """)
    assert result[0] == []
    assert all(value != "accepted" for value in result[1:])


def test_node_visibility_defaults_to_shown_and_requires_a_boolean(tmp_path):
    result = _run(tmp_path, """
        const state = {activeViewName: 'tree'};
        const metadataColumn = name => name === 'group';
        const style = {id: 'node-1', color_by: 'group', shape_by: ''};
        const results = [validateAppearanceStyles([style], 'node')[0].show_node,
            validateAppearanceStyles([{...style, show_node: false}], 'node')[0].show_node];
        try { validateAppearanceStyles([{...style, show_node: 'false'}], 'node'); }
        catch (error) { results.push(error.message); }
        process.stdout.write(JSON.stringify(results));
    """)
    assert result == [True, False, "invalid node visibility"]
