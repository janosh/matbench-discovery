"""Tests for building DFT diatomic reference curves."""

from __future__ import annotations

from typing import Any

import pytest

from matbench_discovery.metrics.diatomics.reference import (
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
    """A one-point deep adiabatic spin drop is replaced by the neighbor spin."""
    candidate_points: CandidatePoints = {
        0: [make_point(1, 0), make_point(2, 0.05), make_point(3, 0)],
        6: [make_point(1, 5), make_point(2, -20, 6), make_point(3, 5)],
    }

    merged, replacements = merge_postprocessed_min_energy_curve(candidate_points)

    assert [point["spin_candidate"] for point in merged] == ["0", "0", "0"]
    assert [point["energies"][-1] for point in merged] == [0, 0.05, 0]
    assert merged[1]["dropped_spin_candidate"] == "6"
    assert replacements[0].kind == "isolated_spin_branch_drop"
    assert replacements[0].distance == 2
    assert replacements[0].original_spin == "6"
    assert replacements[0].original_energy == -20
    assert replacements[0].replacement_energy == 0.05


def test_replaces_drop_between_different_neighbor_spin_branches() -> None:
    """A one-point drop between two smooth neighbor spins uses the smoother branch."""
    candidate_points: CandidatePoints = {
        0: [make_point(1, 0), make_point(2, 1), make_point(3, 5)],
        6: [make_point(1, 5), make_point(2, -20, 6), make_point(3, 5)],
        10: [make_point(1, 5), make_point(2, 0, 10), make_point(3, 0, 10)],
    }

    merged, replacements = merge_postprocessed_min_energy_curve(candidate_points)

    assert [point["spin_candidate"] for point in merged] == ["0", "10", "10"]
    assert [point["energies"][-1] for point in merged] == [0, 0, 0]
    assert replacements[0].kind == "isolated_spin_branch_drop"
    assert replacements[0].original_spin == "6"
    assert replacements[0].replacement_spin == "10"


def test_replaces_low_gain_spin_branch_island() -> None:
    """A tiny isolated adiabatic branch island is replaced by surrounding spin."""
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

    assert [point["spin_candidate"] for point in merged] == ["2", "2", "2", "2"]
    assert [point["energies"][-1] for point in merged] == [0, -0.96, -0.9, 0]
    assert [edit.kind for edit in edits] == [
        "low_gain_spin_branch_island",
        "low_gain_spin_branch_island",
    ]
    assert [edit.expected_energy for edit in edits] == [-0.475, -0.5]


def test_keeps_spin_branch_island_when_replacement_is_less_smooth() -> None:
    """A low-gain spin island is kept if replacement worsens local smoothness."""
    candidate_points: CandidatePoints = {
        2: [make_point(1, 0, 2), make_point(2, -1.05, 2), make_point(3, 0, 2)],
        4: [make_point(1, -1, 4), make_point(2, -0.9, 4), make_point(3, -1, 4)],
    }

    merged, edits = merge_postprocessed_min_energy_curve(candidate_points)

    assert edits == []
    assert [point["spin_candidate"] for point in merged] == ["4", "2", "4"]
    assert [point["energies"][-1] for point in merged] == [-1, -1.05, -1]


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
