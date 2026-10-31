"""Tests for diatomic molecule generation and energy/force calculations."""

import copy

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import all_changes
from ase.calculators.emt import EMT
from ase.formula import Formula

from matbench_discovery.diatomics import (
    DiatomicResults,
    calc_diatomic_curve,
    generate_diatomics,
)


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
        assert atoms.get_chemical_formula() == Formula(formula).format("hill")
        assert atoms.get_distance(0, 1) == pytest.approx(dist)
        # large periodic box so cell-requiring calculators work; images don't interact
        assert all(atoms.pbc)
        assert atoms.cell.lengths() == pytest.approx([50, 50, 50])


def test_generate_diatomics_rejects_distance_beyond_half_box() -> None:
    """Separations >= box_size/2 are rejected (wrong min-image distance under PBC)."""
    with pytest.raises(ValueError, match="must be < box_size/2"):
        generate_diatomics("H", "H", [3.0], box_size=5.0)


class CuRangeLimitedEMT(EMT):
    """EMT raising for Cu pairs beyond cu_max_dist, i.e. a mid-curve failure."""

    def __init__(self, cu_max_dist: float) -> None:
        """Store the Cu separation (Å) above which calculations raise."""
        super().__init__()
        self.cu_max_dist = cu_max_dist

    def calculate(
        self,
        atoms: Atoms | None = None,
        properties: tuple[str, ...] = ("energy",),
        system_changes: list[str] = all_changes,
    ) -> None:
        """Raise RuntimeError for long-range Cu pairs, else defer to EMT."""
        if (
            atoms is not None
            and "Cu" in atoms.get_chemical_symbols()
            and atoms.get_distance(0, 1) > self.cu_max_dist
        ):
            raise RuntimeError(f"synthetic failure beyond {self.cu_max_dist} Å")
        super().calculate(atoms, properties, system_changes)


STALE_CU_CURVE = {
    "distances": [2.0],
    "energies": [-1.0],
    "forces": [[[0.1, 0, 0], [-0.1, 0, 0]]],
}


@pytest.mark.parametrize(
    ("distances", "initial_results", "cu_max_dist", "expected_cu"),
    [
        ([1.0, 2.0], {}, np.inf, "computed"),
        ([1.0], {"Cu-Cu": STALE_CU_CURVE}, np.inf, "computed"),
        # cached samples on the requested grid are kept as is, not recomputed
        ([2.0], {"Cu-Cu": STALE_CU_CURVE}, np.inf, "cached"),
        # a failure after the first point discards that point and any stale cache
        ([1.0, 2.0], {"Cu-Cu": STALE_CU_CURVE}, 1.5, "failed"),
    ],
    ids=["fresh", "stale-distance-grid", "cache-hit", "failed-element-reset"],
)
def test_calc_diatomic_curve_results(
    distances: list[float],
    initial_results: DiatomicResults,
    cu_max_dist: float,
    expected_cu: str,
) -> None:
    """Pairs give full curves, reuse matching caches, replace stale, reset failures."""
    # mixed symbol/number pair specs, incl. H-H requested in both forms
    pairs = [(1, 1), ("H", 1), (8, "O"), ("Cu", 29)]
    initial_results = copy.deepcopy(initial_results)  # results are updated in-place
    calculator = CuRangeLimitedEMT(cu_max_dist)
    results = calc_diatomic_curve(pairs, calculator, "test", distances, initial_results)
    fresh_results = calc_diatomic_curve(pairs, EMT(), "test", distances, {})

    assert set(results) == {"H-H", "O-O", "Cu-Cu"}
    expected_cu_curves = {
        "cached": STALE_CU_CURVE,
        "failed": {"distances": distances, "energies": [], "forces": []},
    }
    assert results["Cu-Cu"] == expected_cu_curves.get(
        expected_cu, fresh_results["Cu-Cu"]
    )
    n_atoms, n_dims = 2, 3
    for formula, curve in fresh_results.items():
        assert curve["distances"] == distances, formula
        assert all(isinstance(energy, int | float) for energy in curve["energies"])
        assert np.array(curve["energies"]).shape == (len(distances),), formula
        expected_shape = (len(distances), n_atoms, n_dims)
        assert np.array(curve["forces"]).shape == expected_shape, formula
        if formula != "Cu-Cu":
            assert results[formula] == curve, formula


def test_calc_diatomic_curve_energy_minimum_and_force_directions() -> None:
    """The energy minimum is interior, forces repel at short and attract at long range
    along the bond axis only, flipping sign at the minimum.
    """
    distances = np.linspace(5, 1.5, 36)  # Å, spans the Cu-Cu equilibrium distance
    results = calc_diatomic_curve([("Cu", "Cu")], EMT(), "test", distances, {})

    energies = results["Cu-Cu"]["energies"]
    min_energy_idx = int(np.argmin(energies))
    assert 0 < min_energy_idx < len(energies) - 1

    forces = np.asarray(results["Cu-Cu"]["forces"])
    # long range (first distance): atoms attract (first atom pulled toward +x)
    assert forces[0, 0, 0] > 0
    assert forces[0, 1, 0] < 0
    # short range (last distance): atoms repel
    assert forces[-1, 0, 0] < 0
    assert forces[-1, 1, 0] > 0
    # the attractive -> repulsive flip happens next to the energy minimum
    first_repulsive_idx = int(np.argmax(forces[:, 0, 0] < 0))
    assert abs(first_repulsive_idx - min_energy_idx) <= 1
    np.testing.assert_array_equal(forces[..., 1:], 0)  # no off-axis components
