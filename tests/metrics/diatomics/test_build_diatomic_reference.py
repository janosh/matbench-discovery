"""Tests for assembling bundled DFT diatomic references."""

from __future__ import annotations

import gzip
import json
from dataclasses import asdict
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from pathlib import Path

from matbench_discovery.metrics.diatomics.reference import CurvePostprocessEdit
from scripts.evals.build_diatomic_reference import (
    build_reference,
    load_candidate_points,
    serializable_curve,
    write_merged_curve,
)

FINITE_FORCES = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]


def write_curve_points(
    tmp_path: Path, points: list[dict[str, object]]
) -> dict[str, object]:
    """Write test curve points and return the shared valid point."""
    curve_dir = tmp_path / "H_pbe_n0"
    curve_dir.mkdir()
    valid_point = {"distance": 1.0, "energies": [-1.0], "forces": FINITE_FORCES}
    (curve_dir / "curve.json").write_text(
        json.dumps({"points": [valid_point, *points]})
    )
    return valid_point


def test_load_candidate_points_skips_malformed_points(tmp_path: Path) -> None:
    """Malformed VASP points are skipped instead of aborting the load."""
    valid_point = write_curve_points(
        tmp_path,
        [
            {"distance": 1.1, "energies": [-0.9], "forces": bad_forces}
            for bad_forces in (
                [None],
                [[]],
                [[0, 0, 0]],
                [[0, 0, 0], [0, 0]],
                [[0, 0, float("nan")], [0, 0, 0]],
                [[0, 0, "x"], [0, 0, 0]],
            )
        ]
        + [{"distance": "bad", "energies": [-0.8], "forces": FINITE_FORCES}],
    )

    loaded_points = load_candidate_points(str(tmp_path), "H", "pbe", 0)

    assert loaded_points == [valid_point]


@pytest.mark.parametrize(
    ("candidate_map", "expected_skipped", "expected_quality"),
    [
        ({}, {}, []),
        # a pair merging to <2 points (H has one valid point) is counted as skipped
        # and logged in the quality report instead of silently vanishing from the
        # bundled reference
        (
            {"H": [0]},
            {"PBE": 1, "r2SCAN": 1},
            [("H", "pbe", True), ("H", "r2scan", True)],
        ),
    ],
    ids=["empty_map", "too_few_points"],
)
def test_build_reference_reports_skipped_pairs(
    tmp_path: Path,
    candidate_map: dict[str, list[int]],
    expected_skipped: dict[str, int],
    expected_quality: list[tuple[str, str, bool]],
) -> None:
    """Empty/too-short inputs yield an empty reference plus skip diagnostics."""
    candidate_map_path = tmp_path / "candidate-map.json"
    out_path = tmp_path / "reference.json.gz"
    merged_dir = tmp_path / "merged"
    candidate_map_path.write_text(json.dumps(candidate_map))
    if candidate_map:
        write_curve_points(tmp_path, [])  # one valid point -> merged curve too short

    summary = build_reference(
        src_dir=str(tmp_path),
        candidate_map_path=str(candidate_map_path),
        out_path=str(out_path),
        merged_dir=str(merged_dir),
        min_drop_ev=3.0,
        postprocess=True,
    )

    with gzip.open(out_path, "rt", encoding="utf-8") as file:
        assert json.load(file) == {"PBE": {}, "r2SCAN": {}}
    assert summary == {
        "merged": {},
        "skipped": expected_skipped,
        "short_candidate_pairs": {},
        "tail_jump_pairs": {},
        "magmom_jump_pairs": {},
        "postprocess_edits": {},
    }
    quality_rows = json.loads((merged_dir / "reference-quality.json").read_text())
    assert [
        (row["symbol"], row["xc"], row["skipped"]) for row in quality_rows
    ] == expected_quality


def test_build_reference_reports_dissociation_tail_jump(tmp_path: Path) -> None:
    """A large final energy step is recorded without modifying the merged curve."""
    candidate_map_path = tmp_path / "candidate-map.json"
    out_path = tmp_path / "reference.json.gz"
    merged_dir = tmp_path / "merged"
    candidate_map_path.write_text(json.dumps({"H": [0]}))
    write_curve_points(
        tmp_path,
        [
            {"distance": 2.0, "energies": [-0.95], "forces": FINITE_FORCES},
            {"distance": 3.0, "energies": [-0.5], "forces": FINITE_FORCES},
        ],
    )

    summary = build_reference(
        src_dir=str(tmp_path),
        candidate_map_path=str(candidate_map_path),
        out_path=str(out_path),
        merged_dir=str(merged_dir),
        min_drop_ev=3.0,
        postprocess=True,
    )

    assert summary["tail_jump_pairs"] == {"PBE": 1}
    assert summary["magmom_jump_pairs"] == {}
    quality_rows = json.loads((merged_dir / "reference-quality.json").read_text())
    pbe_row = next(row for row in quality_rows if row["xc"] == "pbe")
    assert pbe_row["tail_jumps"] == 1


def test_serializable_curve_reports_magmoms_and_spin_candidates() -> None:
    """Per-point magmoms and winning spin candidates survive serialization."""
    points = [
        {
            "distance": 1.0,
            "energies": [-1.0],
            "forces": FINITE_FORCES,
            "magmoms": [1.5464, -1.5464],
            "spin_candidate": "afm",
        },
        # points missing magmoms (e.g. truncated OUTCAR) serialize as null
        {"distance": 2.0, "energies": [-0.5], "forces": FINITE_FORCES},
    ]

    curve = serializable_curve(points)

    assert curve["magmoms"] == [[1.546, -1.546], None]
    assert curve["spin_candidates"] == ["afm", None]
    assert len(curve["distances"]) == len(curve["magmoms"]) == 2


def test_write_merged_curve_serializes_postprocess_edits(tmp_path: Path) -> None:
    """Merged diagnostics include serialized postprocess edit metadata."""
    edit = CurvePostprocessEdit(
        kind="isolated_energy_bump",
        distance=1.2,
        original_spin="2",
        replacement_spin=None,
        original_energy=-0.5,
        replacement_energy=None,
        expected_energy=-0.8,
    )
    point = {"distance": 1.2, "energies": [-0.5], "forces": FINITE_FORCES}

    write_merged_curve(str(tmp_path), "H", "pbe", [point], [edit])

    curve_data = json.loads((tmp_path / "H_pbe" / "curve.json").read_text())
    assert curve_data["postprocess"]["postprocess_edits"] == [asdict(edit)]
