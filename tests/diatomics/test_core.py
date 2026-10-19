"""Tests for diatomic molecule generation and energy/force calculations."""

import numpy as np
import pytest
from ase import Atoms
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


@pytest.mark.parametrize(
    ("distances", "initial_results"),
    [
        ([1.0, 2.0], {}),
        (
            [1.0],
            {
                "Cu-Cu": {
                    "distances": [2.0],
                    "energies": [-1.0],
                    "forces": [[[0.1, 0, 0], [-0.1, 0, 0]]],
                }
            },
        ),
    ],
    ids=["fresh", "stale-distance-grid"],
)
def test_calc_diatomic_curve_results(
    distances: list[float], initial_results: DiatomicResults
) -> None:
    """Symbol/number pairs yield full curves, replacing caches on stale grids."""
    # mixed symbol/number pair specs, incl. H-H requested in both forms
    pairs = [(1, 1), ("H", 1), (8, "O"), ("Cu", 29)]
    results = calc_diatomic_curve(pairs, EMT(), "test", distances, initial_results)

    assert set(results) == {"H-H", "O-O", "Cu-Cu"}
    assert results["Cu-Cu"]["energies"] != [-1.0]
    n_atoms, n_dims = 2, 3
    for formula, curve in results.items():
        assert curve["distances"] == distances, formula
        assert all(isinstance(energy, int | float) for energy in curve["energies"])
        assert np.array(curve["energies"]).shape == (len(distances),), formula
        expected_shape = (len(distances), n_atoms, n_dims)
        assert np.array(curve["forces"]).shape == expected_shape, formula


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
