"""Tests for the run_diatomics CLI, sharding, and merge behavior."""

import gzip
import json
import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest
from ase.data import atomic_numbers, covalent_radii

from tests.utils import import_repo_script


@pytest.fixture(scope="module", name="run_diatomics")
def run_diatomics_fixture() -> ModuleType:
    """models/run_diatomics.py imported once for all tests in this module."""
    return import_repo_script("run_diatomics", "models/run_diatomics.py")


@pytest.mark.parametrize(
    ("argv_tail", "expected_stdout"),
    [
        ("--model emt --min-dist 0", None),
        ("--model emt --min-dist -0.1", None),
        ("--model emt --max-dist 0.1", None),
        ("--model emt --min-dist 2 --max-dist 1", None),
        ("--model emt --max-z 0", None),
        ("--model emt --n-points 1", None),
        ("--model mace-mp-0", "--model mace_mp_0"),
        ("--model matris_10m_oam", "--model matris_10m_oam"),
        ("--model matris_10m_mp", "--model matris_10m_mp"),
        ("--model sevennet_0", "--model sevennet_0"),
        ("--model sevennet_l3i5", "--model sevennet_l3i5"),
        ("--model sevennet_mf_ompa", "--model sevennet_mf_ompa"),
        ("--model sevennet_omni_i12", "--model sevennet_omni_i12"),
        ("--model alphanet_v1_mptrj", "--model alphanet_v1_mptrj"),
        ("--model alphanet_v1_oam", "--model alphanet_v1_oam"),
        # --print-cmd must forward computation params, not just --model (else a copied
        # command silently runs with defaults)
        ("--model mace-mp-0 --min-dist 0.2", "--min-dist 0.2"),
    ],
)
def test_run_diatomics_cli_validation(
    argv_tail: str,
    expected_stdout: str | None,
    run_diatomics: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """run_diatomics validates inputs and resolves shared-runner models."""
    monkeypatch.setattr(
        sys, "argv", ["run_diatomics", "--print-cmd", *argv_tail.split()]
    )

    if expected_stdout is None:
        with pytest.raises(SystemExit):
            run_diatomics.main()
        return
    assert run_diatomics.main() == 0
    assert expected_stdout in capsys.readouterr().out


def test_trim_curve_to_finite(run_diatomics: ModuleType) -> None:
    """Non-finite samples outside the scored window are trimmed, inside it fatal."""
    distances = np.geomspace(0.1, 6.0, 30).tolist()
    finite_curve = {
        "distances": distances,
        "energies": [1.0] * 30,
        "forces": [[[0.0] * 3] * 2] * 30,
    }
    assert run_diatomics.trim_curve_to_finite("Cu-Cu", finite_curve) == finite_curve

    # Cu wall scoring starts at 0.8 * r_cov ~ 1.06 A; NaN below that gets trimmed
    deep_overlap_nan = {**finite_curve, "energies": [np.nan, *[1.0] * 29]}
    trimmed = run_diatomics.trim_curve_to_finite("Cu-Cu", deep_overlap_nan)
    assert trimmed is not None
    assert len(trimmed["energies"]) == len(trimmed["distances"]) == 29
    assert trimmed["distances"] == distances[1:]

    wall_r_min = 0.8 * covalent_radii[atomic_numbers["Cu"]]
    wall_idx = int(np.flatnonzero(np.asarray(distances) >= wall_r_min)[0])
    wall_nan_energies = [1.0] * 30
    wall_nan_energies[wall_idx] = np.nan
    wall_nan = {**finite_curve, "energies": wall_nan_energies}
    assert run_diatomics.trim_curve_to_finite("Cu-Cu", wall_nan) is None

    # non-finite at scored separations (last point = 6 A is inside Cu's window since
    # 3.1 * r_vdw(Cu) > 6 A) drops the whole curve
    in_window_nan = {**finite_curve, "energies": [*[1.0] * 29, np.inf]}
    assert run_diatomics.trim_curve_to_finite("Cu-Cu", in_window_nan) is None

    # H's window ends at 3.1 * r_vdw(H) ~ 3.7 A; a NaN at 6 A is above it (never
    # scored) so only that point is trimmed, not the whole curve
    above_window_nan = {**finite_curve, "energies": [*[1.0] * 29, np.nan]}
    trimmed = run_diatomics.trim_curve_to_finite("H-H", above_window_nan)
    assert trimmed is not None
    assert len(trimmed["energies"]) == len(trimmed["distances"]) == 29
    assert trimmed["distances"] == distances[:-1]

    empty_curve = {"distances": [], "energies": [], "forces": []}
    assert run_diatomics.trim_curve_to_finite("Cu-Cu", empty_curve) is None


def test_run_diatomics_exclusion_reasons(run_diatomics: ModuleType) -> None:
    """Curated + run-discovered exclusion reasons and their use in metric drops."""
    curated_reasons = {"He-He": "exploding errors"}

    assert (
        run_diatomics.get_excluded_formula_reasons("alphanet_v1_oam") == curated_reasons
    )
    # curated reasons win over run-discovered invalid curves (He-He), unknown
    # formulas fall back to the generic invalid-curve reason (X-X), and non-MP
    # formulas (At-At) are never recorded since metrics skip them benchmark-wide
    assert run_diatomics.get_excluded_formula_reasons(
        "alphanet_v1_oam", ["He-He", "At-At", "X-X"]
    ) == curated_reasons | {"X-X": "invalid or unsupported curve"}

    # drop_metric_exclusions removes excluded formula and element-keyed metrics
    metrics = {"He": {"energy_jump": 1.0}, "Cu": {"energy_jump": 2.0}}
    assert run_diatomics.drop_metric_exclusions("alphanet_v1_oam", metrics) == {
        "Cu": {"energy_jump": 2.0}
    }


def write_diatomics_shard(
    shard_dir: Path,
    z_value: int,
    curves: dict[str, dict[str, list[float] | list[list[list[float]]]]],
    run_metadata: dict[str, object] | None = None,
) -> None:
    """Write a minimal run_diatomics shard file for merge-shards tests."""
    shard_dir.mkdir(parents=True, exist_ok=True)
    shard_path = shard_dir / f"Z{z_value:03d}-diatomics.json.gz"
    shard_data = {
        "homo-nuclear": curves,
        "hetero-nuclear": {},
        "distances": [1.0, 2.0],
    }
    if run_metadata is not None:
        shard_data["run_metadata"] = run_metadata
    with gzip.open(shard_path, mode="wt") as file:
        json.dump(shard_data, file)


@pytest.mark.parametrize(
    ("model_key", "second_shard_metadata", "should_raise"),
    [
        (
            "alphanet_v1_oam",
            {"excluded_formula_reasons": {"H-H": "must not leak"}},
            False,
        ),
        ("emt", {"excluded_formula_reasons": {"He-He": "invalid curve"}}, True),
    ],
)
def test_run_diatomics_merge_shards_exclusion_gate(
    model_key: str,
    second_shard_metadata: dict[str, object],
    should_raise: bool,
    run_diatomics: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """merge-shards accepts missing formulas only when YAML-curated."""
    shard_dir = tmp_path / f"{run_diatomics.today}-diatomics-shards"
    curves = {
        "H-H": {
            "energies": [0.0, 1.0],
            "forces": [[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]] * 2,
        }
    }
    write_diatomics_shard(shard_dir, 1, curves)
    write_diatomics_shard(shard_dir, 2, {}, run_metadata=second_shard_metadata)
    cli = f"run_diatomics --model {model_key} --merge-shards --max-z 2 --out-dir"
    monkeypatch.setattr(sys, "argv", [*cli.split(), str(tmp_path)])

    if should_raise:
        with pytest.raises(SystemExit) as exc_info:
            run_diatomics.main()
        assert exc_info.value.code == 2
        stderr = capsys.readouterr().err
        assert "Missing curves in shards" in stderr
        assert "He-He" in stderr
    else:
        assert run_diatomics.main() == 0
        merged_path = tmp_path / f"{run_diatomics.today}-diatomics.json.gz"
        with gzip.open(merged_path, mode="rt") as file:
            merged_data = json.load(file)
        assert merged_data["run_metadata"]["excluded_formula_reasons"] == {
            "He-He": "exploding errors"
        }
