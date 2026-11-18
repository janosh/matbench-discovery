"""Tests for formation and reference energy calculations."""

from collections.abc import Callable
from typing import Any

import pytest
from matscikit import Composition, Structure

from matbench_discovery.energy import (
    Entry,
    calc_energy_from_e_refs,
    get_elemental_ref_entries,
)


@pytest.fixture
def ref_energies() -> dict[str, float]:
    """Return elemental reference energies for Fe-O test entries."""
    return {"Fe": -1.0, "O": -2.0}


@pytest.mark.parametrize("constructor", [Entry, dict])
@pytest.mark.parametrize("verbose", [True, False])
def test_get_elemental_ref_entries(
    constructor: Callable[..., Entry | dict[str, Any]],
    verbose: bool,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Elemental reference entries are identified; verbose prints progress."""
    entries = [
        ("Fe1 O1", -2.5),
        ("Fe1", -1.0),
        (comp1 := "Fe1", energy1 := -2.0),
        (comp2 := "O1", energy2 := -1.0),
        ("O3", -2.0),
    ]
    elemental_ref_entries = get_elemental_ref_entries(
        [constructor(composition=comp, energy=energy) for comp, energy in entries],
        verbose=verbose,
    )
    expected = {"Fe": Entry(comp1, energy1), "O": Entry(comp2, energy2)}
    assert elemental_ref_entries == expected
    stdout, stderr = capsys.readouterr()
    assert ("Sorting 5 entries with 2 dimensions" in stdout) is verbose
    assert ("Finding elemental reference entries" in stderr) is verbose


def test_get_elemental_ref_entries_missing_terminal_entries() -> None:
    """Missing elemental entries should report the missing element."""
    entries = [Entry("Fe1 O1", -2.5), Entry("Fe1", -1.0)]
    with pytest.raises(ValueError, match=r"Missing terminal entries: \{'O'\}"):
        get_elemental_ref_entries(entries, verbose=False)


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
    "input_obj,total_energy,expected",
    [
        ("FeO", -5.0, -1.0),  # formula string
        (Composition("Fe2O3"), -10.0, -0.4),  # composition
        ("Fe4O6", -20.0, -0.4),  # complex composition
        ("Fe", -2.0, -1.0),  # single atom
        ("O2", -6.0, -1.0),  # diatomic
        (Entry(Composition("FeO"), -5.0), None, -1.0),
        ({"composition": "FeO", "energy": -5.0}, None, -1.0),
    ],
)
def test_calc_energy_from_e_refs_various_inputs(
    input_obj: str | Composition | Entry,
    total_energy: float | None,
    expected: float,
    ref_energies: dict[str, float],
) -> None:
    """Test calculation with various input types and compositions."""
    energy = calc_energy_from_e_refs(input_obj, ref_energies, total_energy)
    assert energy == pytest.approx(expected)


def test_calc_energy_from_e_refs_error_cases(
    dummy_struct: Structure,
    ref_energies: dict[str, float],
) -> None:
    """Test error handling."""
    with pytest.raises(ValueError, match="total_energy can't be None"):
        calc_energy_from_e_refs(dummy_struct, ref_energies)

    with pytest.raises(ValueError, match="Missing reference energies"):
        calc_energy_from_e_refs(dummy_struct, {"Fe": -1.0}, total_energy=-5.0)

    with pytest.raises(TypeError, match="Expected Entry, Structure"):
        calc_energy_from_e_refs([1, 2, 3], ref_energies, total_energy=-5.0)

    with pytest.raises(ValueError, match="Entry dict missing required keys"):
        calc_energy_from_e_refs({"energy": -5.0}, ref_energies)
