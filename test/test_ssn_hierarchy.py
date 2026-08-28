"""Direct tests for the merge-event data layer that the split charts render.

`matrix_report` and the SSN viewer both draw stems at `largest_merge` with one bead per
distinct merge size, and a moving-sum line over `merge_impact`. These tests pin the
relationships between those quantities, which the charts rely on but cannot check.
"""
import numpy as np
import pytest

from domainator.data_matrix import DenseDataMatrix, MaxTree
from domainator.ssn_hierarchy import (
    MERGE_IMPACT_MIN_CHILD,
    MERGE_IMPACT_PRODUCT,
    component_size_summary_by_threshold,
    filter_merge_event_rows,
    merge_event_moving_sum,
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
    filtered = filter_merge_event_rows(rows, max_merge_events=2)
    assert [row["threshold_value"] for row in filtered] == [10.0, 20.0]

    unfiltered_sum = merge_event_moving_sum(rows, window_fraction=1.0, grid_points=3)
    filtered_sum = merge_event_moving_sum(filtered, window_fraction=1.0, grid_points=3)

    # The unfiltered series spans the full threshold range and totals every event.
    assert max(unfiltered_sum["x"]) == pytest.approx(30.0)
    assert max(filtered_sum["x"]) == pytest.approx(20.0)
    assert max(unfiltered_sum["y"]) == 16
    assert max(filtered_sum["y"]) == 15
