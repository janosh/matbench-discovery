"""Tests for diatomic molecule generation and energy/force calculations."""

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.emt import EMT

from matbench_discovery.diatomics import (
    atom_num_symbol_map,
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
        assert set(atoms.get_chemical_formula()) == set(formula)
        assert atoms.get_distance(0, 1) == pytest.approx(dist)
        assert not any(atoms.pbc)  # periodic boundary conditions should be False


@pytest.mark.parametrize(
    "z1,z2,expected_formula",
    [
        (1, 1, "H-H"),  # hydrogen by atomic number
        ("H", 1, "H-H"),  # mixed string and number
        (8, "O", "O-O"),  # oxygen
        ("Cu", 29, "Cu-Cu"),  # copper (supported by EMT)
    ],
)
def test_atomic_number_to_symbol_conversion(
    z1: str | int, z2: str | int, expected_formula: str
) -> None:
    """Test conversion between atomic numbers and chemical symbols."""
    distances = [1.0]
    pairs = [(z1, z2)]
    # Use EMT calculator as a simple test calculator
    calculator = EMT()
    results: dict[str, dict[str, list[float | list[list[float]]]]] = {}

    results = calc_diatomic_curve(pairs, calculator, "test", distances, results)

    assert expected_formula in results
    assert "energies" in results[expected_formula]
    assert "forces" in results[expected_formula]
    assert len(results[expected_formula]["energies"]) == len(distances)
    assert len(results[expected_formula]["forces"]) == len(distances)


def test_atom_num_symbol_map() -> None:
    """Test the atomic number to symbol mapping."""
    assert atom_num_symbol_map[1] == "H"
    assert atom_num_symbol_map[8] == "O"
    assert atom_num_symbol_map[26] == "Fe"
    assert len(atom_num_symbol_map) > 100  # should contain most elements


def test_calc_diatomic_curve_results() -> None:
    """Test that calc_diatomic_curve returns correct energy and force arrays."""
    calculator = EMT()
    distances = [1.0, 2.0]
    pairs = [("Cu", "Cu")]  # Use Cu instead of H (better EMT support)
    results: dict[str, dict[str, list[float | list[list[float]]]]] = {}

    results = calc_diatomic_curve(pairs, calculator, "test", distances, results)

    assert "Cu-Cu" in results

    cu2_results = results["Cu-Cu"]
    assert len(cu2_results["energies"]) == len(distances)
    assert len(cu2_results["forces"]) == len(distances)
    assert all(isinstance(energy, int | float) for energy in cu2_results["energies"])
    assert np.array(cu2_results["energies"]).shape == (len(distances),)
    n_atoms, n_dims = len(pairs[0]), 3
    assert np.array(cu2_results["forces"]).shape == (len(distances), n_atoms, n_dims)


def test_calc_diatomic_curve_energy_trend() -> None:
    """Test that energy follows expected trend with distance."""
    calculator = EMT()
    # Test with a range of distances around the equilibrium bond length for Cu2
    distances = np.linspace(2.0, 5.0, 10)  # Cu-Cu has larger equilibrium distance
    pairs = [("Cu", "Cu")]
    results: dict[str, dict[str, list[float | list[list[float]]]]] = {}

    results = calc_diatomic_curve(pairs, calculator, "test", distances, results)

    energies = results["Cu-Cu"]["energies"]
    # Energy should have a minimum at the equilibrium distance
    min_energy_idx = np.argmin(energies)
    assert min_energy_idx > 0  # not at the shortest distance
    assert min_energy_idx < len(energies) - 1  # not at the longest distance


def test_calc_diatomic_curve_force_directions() -> None:
    """Test that forces point in the expected directions."""
    calculator = EMT()
    # Test at very short and very long distances
    distances = [2.0, 5.0]  # Å (adjusted for Cu-Cu)
    pairs = [("Cu", "Cu")]
    results: dict[str, dict[str, list[float | list[list[float]]]]] = {}

    results = calc_diatomic_curve(pairs, calculator, "test", distances, results)

    forces: list[list[list[float]]] = results["Cu-Cu"]["forces"]  # type: ignore[assignment]
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
    """Test that calc_diatomic_curve recalculates requested pairs."""
    calculator = EMT()
    distances = [1.0]
    pairs = [("Cu", "Cu")]  # Only test Cu-Cu
    # Pre-populate results with Cu2
    initial_results = {
        "Cu-Cu": {"energies": [-1.0], "forces": [[[0.1, 0, 0], [-0.1, 0, 0]]]}
    }

    # Calculate new results but with same pair
    results = calc_diatomic_curve(pairs, calculator, "test", distances, initial_results)  # type: ignore[arg-type]

    # Cu2 results should be recalculated since it's in pairs
    assert len(results["Cu-Cu"]["energies"]) == 1
    assert len(results["Cu-Cu"]["forces"]) == 1
    # Values should be calculated by EMT, not taken from initial_results
    assert isinstance(results["Cu-Cu"]["energies"][0], int | float)
    assert isinstance(results["Cu-Cu"]["forces"][0], list)
    assert len(results["Cu-Cu"]["forces"][0]) == 2  # two atoms
    assert len(results["Cu-Cu"]["forces"][0][0]) == 3  # three components
