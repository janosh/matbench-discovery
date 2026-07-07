"""Tests for assembling bundled DFT diatomic references."""

from __future__ import annotations

import gzip
import json
from dataclasses import asdict
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

from matbench_discovery.metrics.diatomics.reference import CurvePostprocessEdit
from scripts.evals.build_diatomic_reference import (
    build_reference,
    load_candidate_points,
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


def test_build_reference_writes_empty_quality_file(tmp_path: Path) -> None:
    """Empty inputs still create the optional merged diagnostics file."""
    candidate_map_path = tmp_path / "candidate-map.json"
    out_path = tmp_path / "reference.json.gz"
    merged_dir = tmp_path / "merged"
    candidate_map_path.write_text("{}")

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
        "short_candidate_pairs": {},
        "postprocess_edits": {},
    }
    assert json.loads((merged_dir / "reference-quality.json").read_text()) == []


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
