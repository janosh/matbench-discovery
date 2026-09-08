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


def test_calc_diatomic_curve_energy_trend() -> None:
    """Test that energy follows expected trend with distance."""
    distances = np.logspace(1, -1, 40)  # Cu-Cu has larger equilibrium distance
    results = calc_diatomic_curve([("Cu", "Cu")], EMT(), "test", distances, {})

    energies = results["Cu-Cu"]["energies"]
    # energy minimum must sit at the equilibrium distance, i.e. strictly interior
    min_energy_idx = np.argmin(energies)
    assert 0 < min_energy_idx < len(energies) - 1


def test_calc_diatomic_curve_force_directions() -> None:
    """Test that forces point in the expected directions."""
    distances = [2.0, 5.0]  # Å, very short and very long for Cu-Cu
    results = calc_diatomic_curve([("Cu", "Cu")], EMT(), "test", distances, {})

    forces: list[list[list[float]]] = results["Cu-Cu"]["forces"]  # ty: ignore[invalid-assignment]
    # At short distance, forces should point away from each other
    assert forces[0][0][0] < 0  # first atom, x component
    assert forces[0][1][0] > 0  # second atom, x component
    # At long distance, forces should point toward each other
    assert forces[1][0][0] > 0  # first atom, x component
    assert forces[1][1][0] < 0  # second atom, x component
    # y and z components should be zero
    assert all(atom_force[1:] == [0, 0] for frame in forces for atom_force in frame)
