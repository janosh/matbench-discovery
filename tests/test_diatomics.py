"""Tests for diatomic molecule generation and energy/force calculations."""

import gzip
import json
import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.emt import EMT
from conftest import import_repo_script

from matbench_discovery.diatomics import (
    DiatomicResults,
    atom_num_symbol_map,
    calc_diatomic_curve,
    generate_diatomics,
)


@pytest.fixture(scope="module", name="run_diatomics")
def run_diatomics_fixture() -> ModuleType:
    """models/run_diatomics.py imported once for all tests in this module."""
    return import_repo_script("run_diatomics", "models/run_diatomics.py")


@pytest.mark.parametrize(
    "elem1,elem2,distances,expected_formulas",
    [
        ("H", "H", [1.0], ["H2"]),  # single distance
        ("O", "O", [1.0, 2.0], ["O2", "O2"]),  # multiple distances
        ("H", "O", [1.0], ["OH"]),  # hetero-nuclear
        ("Cu", "O", [1.5, 2.0], ["CuO", "CuO"]),  # metal (Cu supported by EMT)
    ],
)
def test_generate_diatomics(
    elem1: str, elem2: str, distances: list[float], expected_formulas: list[str]
) -> None:
    """Test diatomic molecule generation with various element pairs and distances."""
    atoms_list = generate_diatomics(elem1, elem2, distances)

    assert len(atoms_list) == len(distances)
    for atoms, dist, formula in zip(
        atoms_list, distances, expected_formulas, strict=True
    ):
        assert isinstance(atoms, Atoms)
        assert set(atoms.get_chemical_formula()) == set(formula)
        assert atoms.get_distance(0, 1) == pytest.approx(dist)
        # large periodic box so cell-requiring calculators work; images don't interact
        assert all(atoms.pbc)
        assert atoms.cell.lengths() == pytest.approx([50, 50, 50])


def test_generate_diatomics_rejects_distance_beyond_half_box() -> None:
    """Separations >= box_size/2 are rejected (wrong min-image distance under PBC)."""
    with pytest.raises(ValueError, match="must be < box_size/2"):
        generate_diatomics("H", "H", [3.0], box_size=5.0)


def test_atom_num_symbol_map() -> None:
    """Test the atomic number to symbol mapping."""
    assert atom_num_symbol_map[1] == "H"
    assert atom_num_symbol_map[8] == "O"
    assert atom_num_symbol_map[26] == "Fe"
    assert len(atom_num_symbol_map) > 100  # should contain most elements


def test_calc_diatomic_curve_results() -> None:
    """Pairs given as symbols and/or atomic numbers yield correctly shaped curves."""
    distances = [1.0, 2.0]
    # mixed symbol/number pair specs, incl. H-H requested in both forms
    pairs = [(1, 1), ("H", 1), (8, "O"), ("Cu", 29)]
    results = calc_diatomic_curve(pairs, EMT(), "test", distances, {})

    assert set(results) == {"H-H", "O-O", "Cu-Cu"}
    n_atoms, n_dims = 2, 3
    for formula, curve in results.items():
        assert curve["distances"] == distances, formula
        assert all(isinstance(energy, int | float) for energy in curve["energies"])
        assert np.array(curve["energies"]).shape == (len(distances),), formula
        expected_shape = (len(distances), n_atoms, n_dims)
        assert np.array(curve["forces"]).shape == expected_shape, formula


def test_calc_diatomic_curve_energy_trend() -> None:
    """Test that energy follows expected trend with distance."""
    # Test with a range of distances around the equilibrium bond length for Cu2
    distances = np.logspace(1, -1, 40)  # Cu-Cu has larger equilibrium distance
    results = calc_diatomic_curve([("Cu", "Cu")], EMT(), "test", distances, {})

    energies = results["Cu-Cu"]["energies"]
    # Energy should have a minimum at the equilibrium distance
    min_energy_idx = np.argmin(energies)
    assert min_energy_idx > 0  # not at the shortest distance
    assert min_energy_idx < len(energies) - 1  # not at the longest distance


def test_calc_diatomic_curve_force_directions() -> None:
    """Test that forces point in the expected directions."""
    # Test at very short and very long distances
    distances = [2.0, 5.0]  # Å (adjusted for Cu-Cu)
    results = calc_diatomic_curve([("Cu", "Cu")], EMT(), "test", distances, {})

    forces: list[list[list[float]]] = results["Cu-Cu"]["forces"]  # ty: ignore[invalid-assignment]
    # At short distance, forces should point away from each other
    assert forces[0][0][0] < 0  # first atom, x component
    assert forces[0][1][0] > 0  # second atom, x component
    # At long distance, forces should point toward each other
    assert forces[1][0][0] > 0  # first atom, x component
    assert forces[1][1][0] < 0  # second atom, x component
    # y and z components should be zero
    assert all(f[0][1:] == [0, 0] for f in forces)  # first atom
    assert all(f[1][1:] == [0, 0] for f in forces)  # second atom


def test_calc_diatomic_curve_prior_results() -> None:
    """calc_diatomic_curve recalculates requested pairs, overwriting stale results."""
    distances = [1.0]
    initial_results: DiatomicResults = {
        "Cu-Cu": {
            "distances": [2.0],
            "energies": [-1.0],
            "forces": [[[0.1, 0, 0], [-0.1, 0, 0]]],
        }
    }

    results = calc_diatomic_curve(
        [("Cu", "Cu")], EMT(), "test", distances, initial_results
    )

    # Cu-Cu was recalculated on the new distance grid, not taken from initial_results
    cu_curve = results["Cu-Cu"]
    assert cu_curve["distances"] == distances
    assert cu_curve["energies"] != [-1.0]
    assert np.array(cu_curve["forces"]).shape == (len(distances), 2, 3)


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
    """run_diatomics validates distances and canonicalizes model refs."""
    argv = ["run_diatomics", "--print-cmd", *argv_tail.split()]
    monkeypatch.setattr(sys, "argv", argv)

    if expected_stdout is None:
        with pytest.raises(SystemExit):
            run_diatomics.main()
    else:
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

    # Cu window starts at 0.9 * r_cov ~ 1.19 A; NaN below that gets trimmed
    deep_overlap_nan = {**finite_curve, "energies": [np.nan, *[1.0] * 29]}
    trimmed = run_diatomics.trim_curve_to_finite("Cu-Cu", deep_overlap_nan)
    assert trimmed is not None
    assert len(trimmed["energies"]) == len(trimmed["distances"]) == 29
    assert trimmed["distances"] == distances[1:]

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
