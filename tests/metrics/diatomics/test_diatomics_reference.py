"""Tests for building DFT diatomic reference curves."""

from __future__ import annotations

from typing import Any

import pytest

from matbench_discovery.metrics.diatomics.reference import (
    drop_collapsed_scf_points,
    merge_postprocessed_min_energy_curve,
)

CurvePoint = dict[str, Any]
CandidatePoints = dict[int | str, list[CurvePoint]]


def make_point(distance: float, energy: float, magmom: float = 0) -> CurvePoint:
    """Build a minimal finite VASP curve point."""
    return {
        "distance": distance,
        "energies": [energy],
        "forces": [[0, 0, 0], [0, 0, 0]],
        "magmoms": [magmom, magmom],
    }


def test_replaces_isolated_spin_branch_energy_drop() -> None:
    """A one-point deep adiabatic spin drop is replaced by the neighbor spin.

    The drop is kept below the 20 eV collapsed-SCF threshold so it exercises the
    isolated-spin-branch-drop path, not the per-candidate collapse filter.
    """
    candidate_points: CandidatePoints = {
        0: [make_point(1, 0), make_point(2, 0.05), make_point(3, 0)],
        6: [make_point(1, 5), make_point(2, -10, 6), make_point(3, 5)],
    }

    merged, replacements = merge_postprocessed_min_energy_curve(candidate_points)

    assert [point["spin_candidate"] for point in merged] == ["0", "0", "0"]
    assert [point["energies"][-1] for point in merged] == [0, 0.05, 0]
    assert merged[1]["dropped_spin_candidate"] == "6"
    assert replacements[0].kind == "isolated_spin_branch_drop"
    assert replacements[0].distance == 2
    assert replacements[0].original_spin == "6"
    assert replacements[0].original_energy == -10
    assert replacements[0].replacement_energy == 0.05


def test_replaces_drop_between_different_neighbor_spin_branches() -> None:
    """A one-point drop between two smooth neighbor spins uses the smoother branch."""
    candidate_points: CandidatePoints = {
        0: [make_point(1, 0), make_point(2, 1), make_point(3, 5)],
        6: [make_point(1, 5), make_point(2, -10, 6), make_point(3, 5)],
        10: [make_point(1, 5), make_point(2, 0, 10), make_point(3, 0, 10)],
    }

    merged, replacements = merge_postprocessed_min_energy_curve(candidate_points)

    assert [point["spin_candidate"] for point in merged] == ["0", "10", "10"]
    assert [point["energies"][-1] for point in merged] == [0, 0, 0]
    assert replacements[0].kind == "isolated_spin_branch_drop"
    assert replacements[0].original_spin == "6"
    assert replacements[0].replacement_spin == "10"


def test_keeps_low_gain_spin_branch_islands() -> None:
    """Near-degenerate branch flicker is kept: it IS the adiabatic minimum.

    An earlier pipeline substituted <=0.2 eV spin islands with the surrounding branch
    for cosmetic smoothness. That step was removed as unphysical hand-curation: the
    per-distance minimum is reported as computed, and the published per-point magmoms
    + spin_candidates make any remaining branch flicker visible instead of hidden.
    """
    candidate_points: CandidatePoints = {
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
    }

    merged, edits = merge_postprocessed_min_energy_curve(candidate_points)

    assert edits == []
    assert [point["spin_candidate"] for point in merged] == ["2", "4", "4", "2"]
    assert [point["energies"][-1] for point in merged] == [0, -1.0, -0.95, 0]


def test_removes_isolated_upward_energy_bump() -> None:
    """A one-point same-branch upward SCF glitch is removed from the curve."""
    candidate_points: CandidatePoints = {
        4: [
            make_point(1, 0, 4),
            make_point(2, -1, 4),
            make_point(3, -0.2, 4),
            make_point(4, -1.1, 4),
            make_point(5, 0, 4),
        ]
    }

    merged, edits = merge_postprocessed_min_energy_curve(candidate_points)

    assert [point["distance"] for point in merged] == [1, 2, 4, 5]
    assert [point["energies"][-1] for point in merged] == [0, -1, -1.1, 0]
    assert edits[0].kind == "isolated_energy_bump"
    assert edits[0].replacement_spin is None


@pytest.mark.parametrize(
    ("candidate_points", "expected_spins", "expected_energies"),
    [
        (
            {
                0: [make_point(1, 0), make_point(2, 1), make_point(3, 5)],
                6: [make_point(1, 5), make_point(2, -20), make_point(3, 0)],
            },
            ["0", "6", "6"],
            [0, -20, 0],
        ),
        (
            {
                0: [make_point(1, 0), make_point(2, -2), make_point(3, 0)],
                6: [make_point(1, 5), make_point(2, -4), make_point(3, 5)],
            },
            ["0", "6", "0"],
            [0, -4, 0],
        ),
    ],
    ids=["neighbor_spins_differ", "drop_below_threshold"],
)
def test_keeps_non_isolated_or_small_spin_branch_drops(
    candidate_points: CandidatePoints,
    expected_spins: list[str],
    expected_energies: list[float],
) -> None:
    """Only isolated severe spin-branch drops are postprocessed."""
    merged, replacements = merge_postprocessed_min_energy_curve(candidate_points)

    assert replacements == []
    assert [point["spin_candidate"] for point in merged] == expected_spins
    assert [point["energies"][-1] for point in merged] == expected_energies


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


def test_merge_excludes_collapsed_points_from_min_merge() -> None:
    """A collapsed candidate point cannot win the per-distance minimum merge."""
    candidate_points: CandidatePoints = {
        0: [make_point(1, -74.0), make_point(2, -74.1), make_point(3, -74.2)],
        "afm": [make_point(1, -73.9), make_point(2, -261.0), make_point(3, -74.1)],
    }

    merged, edits = merge_postprocessed_min_energy_curve(candidate_points)

    assert [point["spin_candidate"] for point in merged] == ["0", "0", "0"]
    assert [point["energies"][-1] for point in merged] == [-74.0, -74.1, -74.2]
    assert [edit.kind for edit in edits] == ["collapsed_scf_point"]
    assert edits[0].original_spin == "afm"
    assert edits[0].original_energy == -261.0
