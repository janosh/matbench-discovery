"""Tests for the run_diatomics CLI, sharding, and merge behavior."""

import gzip
import json
import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest
from ase.data import atomic_numbers, covalent_radii

from matbench_discovery.enums import Model
from tests.utils import import_repo_script, patch_model_metadata

# synthetic curated exclusions, decoupled from any real model YAML's live pins
CURATED_MODEL = Model["mace_mp_0"]
CURATED_REASONS = {"He-He": "synthetic curated reason"}


@pytest.fixture
def curated_model(monkeypatch: pytest.MonkeyPatch) -> str:
    """Key of a model whose diatomics metrics carry only CURATED_REASONS."""
    diatomics_metrics = {"excluded_formula_reasons": CURATED_REASONS}
    # keep the real identity fields (model_key etc.) that Model.from_ref scans
    metadata = {**CURATED_MODEL.metadata, "metrics": {"diatomics": diatomics_metrics}}
    patch_model_metadata(monkeypatch, CURATED_MODEL, metadata)
    return CURATED_MODEL.name


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
        # shared-runner models are resolved by name, no per-model branching
        ("--model sevennet_0", "--model sevennet_0"),
        ("--model alphanet_v1_oam", "--model alphanet_v1_oam"),
        # --print-cmd must forward computation params, not just --model (else a copied
        # command silently runs with defaults)
        ("--model mace-mp-0 --min-dist 0.2", "--min-dist 0.2"),
        ("--model emt --checkpoint /ckpt/model.pt", "--checkpoint /ckpt/model.pt"),
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

    wall_r_min = 0.8 * covalent_radii[atomic_numbers["Cu"]]
    wall_idx = int(np.flatnonzero(np.asarray(distances) >= wall_r_min)[0])
    for formula, bad_idx, value, expected_distances in (
        # Cu wall scoring starts at ~1.06 A; a NaN below that is unscored.
        ("Cu-Cu", 0, np.nan, distances[1:]),
        ("Cu-Cu", wall_idx, np.nan, None),
        # 6 A is inside Cu's window but above H's ~3.7 A upper bound.
        ("Cu-Cu", -1, np.inf, None),
        ("H-H", -1, np.nan, distances[:-1]),
    ):
        energies = [1.0] * 30
        energies[bad_idx] = value
        trimmed = run_diatomics.trim_curve_to_finite(
            formula, {**finite_curve, "energies": energies}
        )
        if expected_distances is None:
            assert trimmed is None, (formula, bad_idx)
        else:
            assert trimmed is not None, (formula, bad_idx)
            assert len(trimmed["energies"]) == len(trimmed["distances"]) == 29
            assert trimmed["distances"] == expected_distances

    empty_curve = {"distances": [], "energies": [], "forces": []}
    assert run_diatomics.trim_curve_to_finite("Cu-Cu", empty_curve) is None


def test_run_diatomics_exclusion_reasons(
    run_diatomics: ModuleType, curated_model: str
) -> None:
    """Curated + run-discovered exclusion reasons and their use in metric drops."""
    assert run_diatomics.get_excluded_formula_reasons(curated_model) == CURATED_REASONS
    # curated reasons win over run-discovered invalid curves (He-He), unknown
    # formulas fall back to the generic invalid-curve reason (X-X), and non-MP
    # formulas (At-At) are never recorded since metrics skip them benchmark-wide
    assert run_diatomics.get_excluded_formula_reasons(
        curated_model, ["He-He", "At-At", "X-X"]
    ) == CURATED_REASONS | {"X-X": "invalid or unsupported curve"}

    # drop_metric_exclusions removes excluded formula and element-keyed metrics
    metrics = {"He": {"energy_jump": 1.0}, "Cu": {"energy_jump": 2.0}}
    assert run_diatomics.drop_metric_exclusions(curated_model, metrics) == {
        "Cu": {"energy_jump": 2.0}
    }

    # debug calculators without a Model entry have no curated exclusions, but a
    # mistyped model key fails instead of silently skipping its curated reasons
    assert run_diatomics.get_excluded_formula_reasons("emt", ["X-X"]) == {
        "X-X": "invalid or unsupported curve"
    }
    with pytest.raises(ValueError, match="not found in Model enum"):
        run_diatomics.get_excluded_formula_reasons(f"{curated_model}m")


def write_diatomics_shard(
    shard_dir: Path,
    z_value: int,
    curves: dict[str, dict[str, list[float] | list[list[list[float]]]]],
    run_metadata: dict[str, object] | None = None,
) -> None:
    """Write a minimal run_diatomics shard (one package env unless overridden)."""
    shard_dir.mkdir(parents=True, exist_ok=True)
    shard_path = shard_dir / f"Z{z_value:03d}-diatomics.json.gz"
    shard_data = {
        "homo-nuclear": curves,
        "hetero-nuclear": {},
        "distances": [1.0, 2.0],
    }
    metadata = {"versions": {"ase": "3.29.0"}, "n_processes": 1}  # one per file
    shard_data["run_metadata"] = metadata | (run_metadata or {})
    with gzip.open(shard_path, mode="wt") as file:
        json.dump(shard_data, file)


@pytest.mark.parametrize(
    ("model_key", "second_shard_metadata", "should_raise", "shard_date", "dry_run"),
    [
        (
            CURATED_MODEL.name,
            {"excluded_formula_reasons": {"H-H": "must not leak"}},
            False,
            None,
            False,
        ),
        # shards written on an earlier date are picked up via resolve_sharded_prefix
        # and the merged artifact inherits their date
        (
            CURATED_MODEL.name,
            {"excluded_formula_reasons": {"H-H": "must not leak"}},
            False,
            "2020-01-02",
            False,
        ),
        # dry runs merge only *-dry-run-shards and write into dry-run/, so a smoke test
        # can neither read nor clobber a real run's shards or artifact
        (
            CURATED_MODEL.name,
            {"excluded_formula_reasons": {"H-H": "must not leak"}},
            False,
            None,
            True,
        ),
        (
            "emt",
            {"excluded_formula_reasons": {"He-He": "invalid curve"}},
            True,
            None,
            False,
        ),
    ],
    ids=["curated", "prior-date-shards", "dry-run-isolated", "uncurated-raises"],
)
@pytest.mark.usefixtures("curated_model")
def test_run_diatomics_merge_shards_exclusion_gate(
    model_key: str,
    second_shard_metadata: dict[str, object],
    should_raise: bool,
    shard_date: str | None,
    dry_run: bool,
    run_diatomics: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """merge-shards accepts missing formulas only when YAML-curated."""
    shard_date = shard_date or run_diatomics.today
    dry_suffix = "-dry-run" if dry_run else ""
    shard_dir = tmp_path / f"{shard_date}-diatomics{dry_suffix}-shards"
    curves = {
        "H-H": {
            "energies": [0.0, 1.0],
            "forces": [[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]] * 2,
        }
    }
    write_diatomics_shard(shard_dir, 1, curves)
    write_diatomics_shard(shard_dir, 2, {}, run_metadata=second_shard_metadata)
    cli = f"run_diatomics --model {model_key} --merge-shards --max-z 2"
    if dry_run:  # --dry-run forces max-z 3, so Li-Li must be present too
        write_diatomics_shard(shard_dir, 3, {"Li-Li": curves["H-H"]})
        # a real run's shards from the same day must be left alone
        write_diatomics_shard(tmp_path / f"{shard_date}-diatomics-shards", 1, curves)
        cli += " --dry-run"
    monkeypatch.setattr(sys, "argv", [*cli.split(), "--out-dir", str(tmp_path)])

    if should_raise:
        with pytest.raises(SystemExit) as exc_info:
            run_diatomics.main()
        assert exc_info.value.code == 2
        stderr = capsys.readouterr().err
        assert "Missing curves in shards" in stderr
        assert "He-He" in stderr
    else:
        assert run_diatomics.main() == 0
        merged_name = f"{shard_date}-diatomics.json.gz"
        merged_path = tmp_path / ("dry-run" if dry_run else "") / merged_name
        if dry_run:
            assert not (tmp_path / merged_name).is_file()
        with gzip.open(merged_path, mode="rt") as file:
            merged_data = json.load(file)
        merged_metadata = merged_data["run_metadata"]
        assert merged_metadata["excluded_formula_reasons"] == CURATED_REASONS
        assert merged_metadata["n_processes"] == (3 if dry_run else 2)


@pytest.mark.parametrize("shard_args", [(), ("--n-shards", "2", "--shard-index", "0")])
def test_run_diatomics_fails_without_curves(
    shard_args: tuple[str, ...],
    run_diatomics: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """All-failed runs exit non-zero naming the model and write no artifact."""

    def all_pairs_fail(
        _calculator: object, _model_key: str, z_group: list[int], _distances: object
    ) -> tuple[dict[str, object], list[str]]:
        """Stand-in for predict_curves where every pair raised."""
        return {}, [f"Z{z_value}" for z_value in z_group]

    monkeypatch.setattr(run_diatomics, "predict_curves", all_pairs_fail)
    argv = ["run_diatomics", "--model", "emt", "--dry-run", "--out-dir", str(tmp_path)]
    monkeypatch.setattr(sys, "argv", [*argv, *shard_args])
    with pytest.raises(RuntimeError, match="emt wrote 0"):
        run_diatomics.main()
    assert not list((tmp_path / "dry-run").glob("*-diatomics.json.gz"))
