"""Direct tests for the merge-event data layer that the split charts render.

`matrix_report` and the Domainator Similarity Network Viewer both draw stems at
`largest_merge` with one bead per distinct merge size, and a moving-sum line over
`merge_impact`. These tests pin the relationships between those quantities, which
the charts rely on but cannot check.
"""
import numpy as np
import pytest

from domainator.data_matrix import DenseDataMatrix, MaxTree
from domainator.ssn_hierarchy import (
    COMPONENT_LARGEST_COL,
    MERGE_EVENT_DENSITY_BINS,
    MERGE_IMPACT_MIN_CHILD,
    MERGE_IMPACT_PRODUCT,
    component_size_summary_by_threshold,
    filter_merge_event_rows,
    merge_event_density_bin,
    merge_event_moving_sum,
    merge_event_rank_key,
    threshold_merge_event_rows,
)


def _tie_group_tree():
    """Six nodes whose three strongest edges all score 10.0.

    That single tie group therefore contains three *separate* one-node merges, which is
    exactly the case the split chart redesign exists to distinguish: `merge_impact` is 3
    while the largest individual merge is only 1.
    """
    labels = [f"n{i}" for i in range(6)]
    data = np.zeros((6, 6))
    for i, j, value in [(0, 1, 10.0), (2, 3, 10.0), (4, 5, 10.0), (1, 2, 5.0), (3, 4, 3.0)]:
        data[i, j] = value
        data[j, i] = value
    return MaxTree(DenseDataMatrix(data, labels, labels, data_type="score"))


def _event_rows(metric=MERGE_IMPACT_MIN_CHILD, tree=None):
    tree = _tie_group_tree() if tree is None else tree
    return threshold_merge_event_rows(
        component_size_summary_by_threshold(tree, merge_impact_metric=metric)
    )


def test_component_summary_tracks_running_largest_component():
    summary = component_size_summary_by_threshold(_tie_group_tree())

    assert summary[:, COMPONENT_LARGEST_COL].tolist() == [1.0, 2.0, 2.0, 2.0, 4.0, 6.0]


@pytest.mark.parametrize("metric", [MERGE_IMPACT_MIN_CHILD, MERGE_IMPACT_PRODUCT])
def test_merge_size_counts_invariants(metric):
    rows = _event_rows(metric)
    assert len(rows) > 0

    for row in rows:
        counts = row["merge_size_counts"]
        assert all(isinstance(size, int) for size in counts), "sizes must key as int"
        assert all(count > 0 for count in counts.values())

        assert sum(counts.values()) == row["merge_count"]
        assert max(counts, default=0) == row["largest_merge"]
        # Exact, by construction: both come from the same summary rows.
        assert sum(size * count for size, count in counts.items()) == row["merge_impact"]
        assert row["merge_count"] <= row["summary_row_to"] - row["summary_row_from"] + 1
        assert row["largest_merge"] <= row["merge_impact"]
        if row["merge_count"] == 1:
            assert row["largest_merge"] == row["merge_impact"]


def test_tie_group_separates_largest_merge_from_the_sum():
    """The regression the redesign is for: plotting the sum hides the split sizes."""
    rows = {row["threshold_value"]: row for row in _event_rows()}

    tie_group = rows[10.0]
    assert tie_group["merge_size_counts"] == {1: 3}
    assert tie_group["merge_count"] == 3
    assert tie_group["merge_impact"] == 3.0
    # A stem drawn at merge_impact would be 3x too tall for the clusters it represents.
    assert tie_group["largest_merge"] == 1.0

    for threshold in (5.0, 3.0):
        assert rows[threshold]["merge_size_counts"] == {2: 1}
        assert rows[threshold]["largest_merge"] == rows[threshold]["merge_impact"] == 2.0


def test_zero_impact_rows_are_not_counted_as_merges():
    """A summary row with zero impact joined nothing and must not become a bead.

    `component_size_summary_by_threshold` writes 0.0 whenever an edge's endpoints were
    already in one component. A spanning tree never contains such an edge, so this is
    driven straight through `threshold_merge_event_rows` rather than through a MaxTree.
    """
    # Columns: threshold, largest, avg_non_singleton, merge_impact, and the two deltas.
    summary = np.zeros((4, 6), dtype=float)
    summary[0] = [float("inf"), 0, 0, 0, 0, 0]
    summary[1] = [7.0, 2, 2, 2.0, 0, 0]
    summary[2] = [7.0, 2, 2, 0.0, 0, 0]  # joined nothing
    summary[3] = [7.0, 5, 5, 3.0, 0, 0]

    (row,) = threshold_merge_event_rows(summary)

    assert row["merge_size_counts"] == {2: 1, 3: 1}
    assert row["merge_count"] == 2
    assert row["largest_merge"] == 3.0
    assert row["merge_impact"] == 5.0
    # Three summary rows collapsed into two counted merges.
    assert row["merge_count"] < row["summary_row_to"] - row["summary_row_from"] + 1


def test_non_integral_merge_sizes_do_not_collide():
    """Keying the histogram with int() would truncate 2.4 and 2.6 onto the same bead."""
    summary = np.zeros((3, 6), dtype=float)
    summary[0] = [float("inf"), 0, 0, 0, 0, 0]
    summary[1] = [7.0, 2, 2, 2.4, 0, 0]
    summary[2] = [7.0, 5, 5, 2.6, 0, 0]

    (row,) = threshold_merge_event_rows(summary)

    assert sorted(row["merge_size_counts"]) == pytest.approx([2.4, 2.6])
    assert row["merge_count"] == 2
    assert row["largest_merge"] == pytest.approx(2.6)


def test_moving_sum_degenerate_inputs():
    assert merge_event_moving_sum([]) == {"window": 0.0, "x": [], "y": []}
    # One threshold leaves no range to slide a window over.
    single = [{"threshold_value": 5.0, "merge_impact": 3.0}]
    assert merge_event_moving_sum(single) == {"window": 0.0, "x": [], "y": []}
    infinite = [{"threshold_value": float("inf"), "merge_impact": 3.0}]
    assert merge_event_moving_sum(infinite) == {"window": 0.0, "x": [], "y": []}


def test_moving_sum_is_a_centred_window_over_merge_impact():
    rows = [
        {"threshold_value": 10.0, "merge_impact": 3.0},
        {"threshold_value": 10.5, "merge_impact": 1.0},
        {"threshold_value": 20.0, "merge_impact": 5.0},
        {"threshold_value": 30.0, "merge_impact": 2.0},
    ]
    result = merge_event_moving_sum(rows, window_fraction=0.5, grid_points=5)

    assert result["window"] == pytest.approx(10.0)  # 0.5 * (30 - 10)
    assert result["x"] == pytest.approx([10.0, 15.0, 20.0, 25.0, 30.0])
    # Window is inclusive at both ends: at g=15 it spans [10, 20] and catches all three.
    assert result["y"] == [4, 9, 5, 7, 2]


def test_moving_sum_must_be_computed_before_filtering():
    """`filter_merge_event_rows` ranks by merge_impact, so it drops precisely the small
    events the moving-sum line exists to reveal. Computing the sum afterwards undercounts."""
    rows = [
        {"threshold_value": 10.0, "merge_impact": 8.0, "delta_largest": 0.0,
         "delta_avg_non_singleton": 0.0, "edge_index": 0},
        {"threshold_value": 20.0, "merge_impact": 7.0, "delta_largest": 0.0,
         "delta_avg_non_singleton": 0.0, "edge_index": 1},
        {"threshold_value": 30.0, "merge_impact": 1.0, "delta_largest": 0.0,
         "delta_avg_non_singleton": 0.0, "edge_index": 2},
    ]
    # density_bins=0 isolates the cap: the back-fill would otherwise restore t=30.0
    # precisely because it is the only event in its band, which is a different
    # property (see the density tests below).
    filtered = filter_merge_event_rows(rows, max_merge_events=2, density_bins=0)
    assert [row["threshold_value"] for row in filtered] == [10.0, 20.0]

    unfiltered_sum = merge_event_moving_sum(rows, window_fraction=1.0, grid_points=3)
    filtered_sum = merge_event_moving_sum(filtered, window_fraction=1.0, grid_points=3)

    # The unfiltered series spans the full threshold range and totals every event.
    assert max(unfiltered_sum["x"]) == pytest.approx(30.0)
    assert max(filtered_sum["x"]) == pytest.approx(20.0)
    assert max(unfiltered_sum["y"]) == 16
    assert max(filtered_sum["y"]) == 15


# ---------------------------------------------------------------------------
# Density back-fill: the cap must not leave stretches of the axis blank
# ---------------------------------------------------------------------------


def _event_row(edge_index, threshold_value, merge_impact):
    return {
        "edge_index": edge_index,
        "threshold_value": threshold_value,
        "merge_impact": merge_impact,
        "delta_largest": 0.0,
        "delta_avg_non_singleton": 0.0,
    }


def _crowded_rows(count=400):
    """Big events crammed into the top of the axis, tiny ones spread below it.

    This is the shape that motivated the back-fill: on a connected MST-kNN graph the
    weak tail is single outliers attaching to the giant component, so ranking by
    impact alone keeps only the right-hand end of the axis.
    """
    rows = []
    for index in range(count):
        # Thresholds 0.95..1.00 -- all inside the top 5% band.
        rows.append(_event_row(index, 0.95 + (0.05 * index / count), 100.0 + index))
    for index in range(count, count + 100):
        # Thresholds 0.00..0.95, impact 1: every one of these loses to every row above.
        rows.append(_event_row(index, 0.95 * (index - count) / 100.0, 1.0))
    return rows


def test_density_backfill_covers_every_band_the_cap_emptied():
    rows = _crowded_rows()
    filtered = filter_merge_event_rows(rows, max_merge_events=10)

    lo = min(row["threshold_value"] for row in rows)
    hi = max(row["threshold_value"] for row in rows)
    covered = {
        merge_event_density_bin(row["threshold_value"], lo, hi, MERGE_EVENT_DENSITY_BINS)
        for row in filtered
    }
    # Every band holds at least one of the 500 input rows, so every band must end up
    # represented -- the whole point of the pass.
    assert covered == set(range(MERGE_EVENT_DENSITY_BINS))
    # Without it the cap would have kept only the top band.
    capped_only = filter_merge_event_rows(rows, max_merge_events=10, density_bins=0)
    assert {
        merge_event_density_bin(row["threshold_value"], lo, hi, MERGE_EVENT_DENSITY_BINS)
        for row in capped_only
    } == {MERGE_EVENT_DENSITY_BINS - 1}


def test_density_backfill_is_bounded_by_the_band_count():
    """Between N and N + bins, which is what makes the output size predictable."""
    rows = _crowded_rows()
    for cap in (1, 10, 100):
        filtered = filter_merge_event_rows(rows, max_merge_events=cap)
        assert cap <= len(filtered) <= cap + MERGE_EVENT_DENSITY_BINS


def test_density_backfill_never_displaces_a_top_ranked_event():
    """A back-filled row is an addition, never a replacement."""
    rows = _crowded_rows()
    cap = 10
    expected_top = sorted(rows, key=merge_event_rank_key)[:cap]
    filtered = filter_merge_event_rows(rows, max_merge_events=cap)
    kept = {row["edge_index"] for row in filtered}
    assert {row["edge_index"] for row in expected_top} <= kept


def test_density_backfill_picks_the_strongest_event_in_an_empty_band():
    rows = [
        _event_row(0, 1.00, 500.0),   # top band, kept by the cap
        _event_row(1, 0.10, 3.0),     # low band: three candidates, this is the strongest
        _event_row(2, 0.11, 2.0),
        _event_row(3, 0.12, 1.0),
        _event_row(4, 0.00, 4.0),     # a different low band
    ]
    filtered = filter_merge_event_rows(rows, max_merge_events=1)
    kept = {row["edge_index"] for row in filtered}
    assert 0 in kept and 4 in kept
    # Rows 1-3 share a band; only the strongest of them is added.
    assert kept & {1, 2, 3} == {1}


def test_filter_returns_everything_when_nothing_was_dropped():
    """No cap, or fewer rows than the cap: no back-fill either, nothing was lost."""
    rows = [_event_row(index, index / 10.0, 1.0) for index in range(5)]
    assert filter_merge_event_rows(rows, max_merge_events=0) == rows
    assert filter_merge_event_rows(rows, max_merge_events=5) == rows
    assert filter_merge_event_rows(rows, max_merge_events=50) == rows


def test_filtered_rows_stay_in_edge_index_order():
    """The series is plotted and turned into slider stops in axis order."""
    filtered = filter_merge_event_rows(_crowded_rows(), max_merge_events=10)
    assert [row["edge_index"] for row in filtered] == sorted(
        row["edge_index"] for row in filtered)


def test_density_bin_places_the_axis_maximum_in_the_last_band():
    """Otherwise `hi` would index one past the end and be reported as its own band."""
    assert merge_event_density_bin(1.0, 0.0, 1.0, 20) == 19
    assert merge_event_density_bin(0.0, 0.0, 1.0, 20) == 0
    assert merge_event_density_bin(0.5, 0.0, 1.0, 20) == 10
    # A degenerate axis has no bands to speak of.
    assert merge_event_density_bin(1.0, 1.0, 1.0, 20) == -1
    assert merge_event_density_bin(0.5, 0.0, 1.0, 0) == -1


def test_density_backfill_is_a_no_op_on_a_single_threshold():
    """A tie group is one point, not an axis, so there are no bands to fill.

    Every event sharing one threshold is the shape `many_cat`-style networks take:
    the cap applies and the back-fill has nothing to add.
    """
    rows = [_event_row(index, 7.0, float(index)) for index in range(50)]
    filtered = filter_merge_event_rows(rows, max_merge_events=3)
    assert len(filtered) == 3
    # The three strongest, returned in edge_index order like every other result.
    assert [row["merge_impact"] for row in filtered] == [47.0, 48.0, 49.0]


def test_filter_rejects_a_negative_band_count():
    with pytest.raises(ValueError):
        filter_merge_event_rows(_crowded_rows(), max_merge_events=10, density_bins=-1)


def test_axis_labels_never_call_a_product_impact_a_node_count():
    """A product impact is a product of two sizes, not a count of anything.

    Both split charts read these strings, so getting it wrong here mislabels the
    axis titles and both hover readouts at once.
    """
    from domainator.ssn_hierarchy import merge_impact_axis_labels

    min_child = merge_impact_axis_labels(MERGE_IMPACT_MIN_CHILD)
    product = merge_impact_axis_labels(MERGE_IMPACT_PRODUCT)

    assert min_child["impact_amount"] == "{} nodes"
    assert min_child["moving_sum_hover"] == "Nodes displaced within"
    assert product["impact_amount"] == "impact {}"
    assert product["moving_sum_hover"] == "Split impact within"

    # Exactly one placeholder, so a caller can fill it with its own interpolation
    # syntax -- Plotly's %{y} in matrix_report, a JS string in the viewer.
    for labels in (min_child, product):
        assert labels["impact_amount"].count("{}") == 1
        assert labels["impact_amount"].format(7) in ("7 nodes", "impact 7")

    # And nothing on the product side says "node".
    assert not any("node" in value.lower() for value in product.values())

    with pytest.raises(ValueError):
        merge_impact_axis_labels("not_a_metric")
