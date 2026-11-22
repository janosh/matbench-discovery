"""Tests for formation and reference energy calculations."""

import numpy as np
import pytest
from matscikit import Composition

from matbench_discovery.energy import Entry, formation_energies


@pytest.fixture
def ref_energies() -> dict[str, float]:
    """Return elemental reference energies for Fe-O test entries."""
    return {"Fe": -1.0, "O": -2.0}


def test_entry_from_dict_applies_corrections() -> None:
    """Serialized entries add their total correction like ComputedEntry.energy."""
    base = {"composition": {"Fe": 2.0}, "energy": -4.0}
    assert Entry.from_dict(base).energy == -4.0
    assert Entry.from_dict(base).energy_per_atom == -2.0
    # ComputedEntry("Fe2O3", -10, energy_adjustments=[CompositionEnergyAdjustment(
    # -0.5, 3), ConstantEnergyAdjustment(0.2)]).as_dict() (pymatgen 2026.9.24):
    # adjustments carry adj_per_atom/n_atoms or value; correction is their sum
    computed_entry = {
        "composition": {"Fe": 2.0, "O": 3.0},
        "energy": -10,
        "correction": -1.3,
        "energy_adjustments": [
            {"@class": "CompositionEnergyAdjustment", "adj_per_atom": -0.5},
            {"@class": "ConstantEnergyAdjustment", "value": 0.2},
        ],
        "entry_id": "mp-19770",
    }
    entry = Entry.from_dict(computed_entry)
    assert (entry.energy, entry.entry_id) == (-11.3, "mp-19770")
    assert entry.composition == Composition("Fe2O3")


@pytest.mark.parametrize(
    "composition,total_energy,expected",
    [
        ("FeO", -5.0, -1.0),  # formula string
        (Composition("Fe2O3"), -10.0, -0.4),
        ({"Fe": 4, "O": 6}, -20.0, -0.4),  # mapping, non-reduced
        ("Fe2 O3", -10.0, -0.4),  # pymatgen-style formula with spaces
        ("Fe", -2.0, -1.0),  # single atom
        ("O2", -6.0, -1.0),  # diatomic
    ],
)
def test_formation_energies(
    composition: str | Composition | dict[str, int],
    total_energy: float,
    expected: float,
    ref_energies: dict[str, float],
) -> None:
    """Total energies become per-atom energies relative to the elemental references."""
    e_forms = formation_energies([composition] * 2, [total_energy] * 2, ref_energies)
    np.testing.assert_allclose(e_forms, [expected] * 2, rtol=0, atol=1e-12)


def test_formation_energies_errors(ref_energies: dict[str, float]) -> None:
    """Elements without a reference energy raise; empty input gives an empty array."""
    with pytest.raises(ValueError, match="no elemental entry for N"):
        formation_energies(["FeN"], [-5.0], ref_energies)
    assert formation_energies([], [], ref_energies).shape == (0,)


@pytest.mark.parametrize("energies", [[-2.0], [[-2.0], [-3.0]], [-2.0, -3.0, -4.0]])
def test_formation_energy_shape(energies: list, ref_energies: dict[str, float]) -> None:
    """Reject a missing, nested or surplus energy before NumPy can broadcast it."""
    with pytest.raises(ValueError, match="one total energy for each of 2 compositions"):
        formation_energies(["Fe", "O"], np.asarray(energies), ref_energies)
