"""Tests for building DFT diatomic reference curves."""

from __future__ import annotations

from typing import Any

import pytest

from matbench_discovery.metrics.diatomics.reference import (
    count_dissociation_tail_jumps,
    count_magmom_discontinuities,
    drop_collapsed_scf_points,
    merge_postprocessed_min_energy_curve,
)

CurvePoint = dict[str, Any]
CandidatePoints = dict[int | str, list[CurvePoint]]


def make_point(distance: float, energy: float, magmom: float = 0) -> CurvePoint:
    """Build a minimal finite DFT curve point."""
    return {
        "distance": distance,
        "energies": [energy],
        "forces": [[0, 0, 0], [0, 0, 0]],
        "magmoms": [magmom, magmom],
    }


@pytest.mark.parametrize(
    ("energies", "expected"),
    [
        ([], 0),
        ([0], 0),
        ([0, 0.02, 0.05], 0),
        ([0, 0.05, 0.3], 1),
        ([0, 0.3, 0], 2),
        ([100, 0, 0.01, 0.02], 0),  # short-range repulsion is not inspected
    ],
)
def test_counts_dissociation_tail_jumps(energies: list[float], expected: int) -> None:
    """Large steps among up to three tail points are counted, not the repulsive wall."""
    assert count_dissociation_tail_jumps(energies) == expected


@pytest.mark.parametrize(
    ("magmoms", "expected"),
    [
        ([[1.0, -1.0], [1.1, -1.1], [1.2, -1.2]], 0),
        ([[1.0, -1.0], [2.2, -1.1], [2.3, -2.4]], 2),
        ([[1.0, -1.0], None, [3.0, -3.0]], 0),
    ],
)
def test_counts_magmom_discontinuities(
    magmoms: list[list[float] | None], expected: int
) -> None:
    """Large per-site moment changes are counted while missing points are skipped."""
    assert count_magmom_discontinuities(magmoms) == expected


@pytest.mark.parametrize(
    ("candidate_points", "expected_points", "expected_edits"),
    [
        # one-point deep adiabatic spin drop is replaced by the neighbor spin; the drop
        # stays below the 20 eV collapsed-SCF threshold so it exercises this path, not
        # the per-candidate collapse filter
        pytest.param(
            {
                0: [make_point(1, 0), make_point(2, 0.05), make_point(3, 0)],
                6: [make_point(1, 5), make_point(2, -10, 6), make_point(3, 5)],
            },
            [(1, "0", 0), (2, "0", 0.05), (3, "0", 0)],
            [("isolated_spin_branch_drop", 2, "6", "0", -10, 0.05)],
            id="isolated_drop",
        ),
        # one-point drop between two smooth neighbor spins uses the smoother branch
        pytest.param(
            {
                0: [make_point(1, 0), make_point(2, 1), make_point(3, 5)],
                6: [make_point(1, 5), make_point(2, -10, 6), make_point(3, 5)],
                10: [make_point(1, 5), make_point(2, 0, 10), make_point(3, 0, 10)],
            },
            [(1, "0", 0), (2, "10", 0), (3, "10", 0)],
            [("isolated_spin_branch_drop", 2, "6", "10", -10, 0)],
            id="drop_between_neighbor_spins",
        ),
        pytest.param(
            {
                0: [make_point(1, 0), make_point(2, 1), make_point(3, 5)],
                6: [make_point(1, 5), make_point(2, -20), make_point(3, 0)],
            },
            [(1, "0", 0), (2, "6", -20), (3, "6", 0)],
            [],
            id="neighbor_spins_differ",  # not isolated, so kept
        ),
        pytest.param(
            {
                0: [make_point(1, 0), make_point(2, -2), make_point(3, 0)],
                6: [make_point(1, 5), make_point(2, -4), make_point(3, 5)],
            },
            [(1, "0", 0), (2, "6", -4), (3, "0", 0)],
            [],
            id="drop_below_threshold",
        ),
        # near-degenerate (<=0.2 eV) branch flicker IS the adiabatic minimum: the
        # per-distance minimum is reported as computed, the published per-point
        # magmoms + spin_candidates make the flicker visible
        pytest.param(
            {
                2: [
                    make_point(1, 0),
                    make_point(2, -0.96, 2),
                    make_point(3, -0.90, 2),
                    make_point(4, 0),
                ],
                4: [
                    make_point(1, 1, 4),
                    make_point(2, -1.00, 4),
                    make_point(3, -0.95, 4),
                    make_point(4, 1, 4),
                ],
            },
            [(1, "2", 0), (2, "4", -1.0), (3, "4", -0.95), (4, "2", 0)],
            [],
            id="low_gain_islands",
        ),
        # one-point same-branch upward SCF glitch is removed from the curve
        pytest.param(
            {
                4: [
                    make_point(dist, energy, 4)
                    for dist, energy in enumerate((0, -1, -0.2, -1.1, 0), start=1)
                ]
            },
            [(1, "4", 0), (2, "4", -1), (4, "4", -1.1), (5, "4", 0)],
            [("isolated_energy_bump", 3, "4", None, -0.2, None)],
            id="upward_bump",
        ),
        # a collapsed candidate point cannot win the per-distance minimum merge
        pytest.param(
            {
                0: [make_point(1, -74.0), make_point(2, -74.1), make_point(3, -74.2)],
                "afm": [
                    make_point(1, -73.9),
                    make_point(2, -261),
                    make_point(3, -74.1),
                ],
            },
            [(1, "0", -74.0), (2, "0", -74.1), (3, "0", -74.2)],
            [("collapsed_scf_point", 2, "afm", None, -261, None)],
            id="collapsed_point",
        ),
    ],
)
def test_merge_postprocessed_min_energy_curve(
    candidate_points: CandidatePoints,
    expected_points: list[tuple[float, str, float]],
    expected_edits: list[tuple[object, ...]],
) -> None:
    """Only isolated severe spin drops, bumps and collapsed points are postprocessed."""
    merged, edits = merge_postprocessed_min_energy_curve(candidate_points)
    by_distance = {point["distance"]: point for point in merged}
    assert [
        (dist, point["spin_candidate"], point["energies"][-1])
        for dist, point in by_distance.items()
    ] == expected_points
    assert [
        (
            edit.kind,
            edit.distance,
            edit.original_spin,
            edit.replacement_spin,
            edit.original_energy,
            edit.replacement_energy,
        )
        for edit in edits
    ] == expected_edits
    for edit in edits:  # replaced spin-branch points record the candidate they dropped
        if edit.kind == "isolated_spin_branch_drop":
            dropped = by_distance[edit.distance]["dropped_spin_candidate"]
            assert dropped == edit.original_spin


@pytest.mark.parametrize(
    ("energies", "expected_kept", "n_dropped"),
    [
        # single collapsed interior point (Gd/r2SCAN-style -261 eV between -74 eV)
        ([-74, -261, -74.2, -74.3], [-74, -74.2, -74.3], 1),
        # two consecutive collapsed points peel off iteratively
        ([-74, -500, -700, -74.2], [-74, -74.2], 2),
        # collapsed endpoint (only one neighbor)
        ([-160, -12, -12.1, -12.2], [-12, -12.1, -12.2], 1),
        # eV-scale wells stay untouched
        ([0, -5, -8, -5, 0], [0, -5, -8, -5, 0], 0),
    ],
    ids=["interior", "consecutive", "endpoint", "physical_well"],
)
def test_drop_collapsed_scf_points(
    energies: list[float], expected_kept: list[float], n_dropped: int
) -> None:
    """Variationally collapsed SCF points are dropped per candidate branch."""
    points = [make_point(idx + 1.0, energy) for idx, energy in enumerate(energies)]

    kept, edits = drop_collapsed_scf_points(points, spin_candidate="8")

    assert [point["energies"][-1] for point in kept] == expected_kept
    assert len(edits) == n_dropped
    assert all(edit.kind == "collapsed_scf_point" for edit in edits)
    assert all(edit.original_spin == "8" for edit in edits)
