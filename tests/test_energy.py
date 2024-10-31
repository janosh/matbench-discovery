from collections.abc import Callable
from typing import Any

import pytest
from pymatgen.analysis.phase_diagram import PDEntry
from pymatgen.core import Composition, Structure
from pymatgen.entries.computed_entries import ComputedEntry, Entry
from pymatgen.util.typing import EntryLike

from matbench_discovery.energy import (
    calc_energy_from_e_refs,
    get_e_form_per_atom,
    get_elemental_ref_entries,
    mp_elem_ref_entries,
    mp_elemental_ref_energies,
)


@pytest.fixture
def ref_energies() -> dict[str, float]:
    return {"Fe": -1.0, "O": -2.0}


def test_get_e_form_per_atom() -> None:
    """Test formation energy calculation."""
    entry = {"composition": {"Fe": 1, "O": 1}, "energy": -2.5}
    ref_energies = {"Fe": -1.0, "O": -1.0}
    assert get_e_form_per_atom(entry, ref_energies) == -0.25


@pytest.mark.parametrize("constructor", [PDEntry, ComputedEntry, lambda **x: x])
@pytest.mark.parametrize("verbose", [True, False])
def test_get_elemental_ref_entries(
    constructor: Callable[..., Entry | dict[str, Any]], verbose: bool
) -> None:
    """Test that elemental reference entries are correctly identified."""
    entries = [
        ("Fe1 O1", -2.5),
        ("Fe1", -1.0),
        ("Fe1", -2.0),
        ("O1", -1.0),
        ("O3", -2.0),
    ]
    elemental_ref_entries = get_elemental_ref_entries(
        [constructor(composition=comp, energy=energy) for comp, energy in entries],
        verbose=verbose,
    )
    if constructor.__name__ == "<lambda>":
        expected = {"Fe": PDEntry(*entries[2]), "O": PDEntry(*entries[3])}
    else:
        expected = {"Fe": constructor(*entries[2]), "O": constructor(*entries[3])}
    assert elemental_ref_entries == expected


def test_mp_ref_energies() -> None:
    """Test MP elemental reference energies are in sync with PDEntries saved to disk."""
    for key, val in mp_elemental_ref_energies.items():
        actual = mp_elem_ref_entries[key].energy_per_atom
        assert actual == pytest.approx(val, abs=1e-3), f"{key=}"


@pytest.mark.parametrize(
    "input_obj,total_energy,expected",
    [
        ("FeO", -5.0, -1.0),  # formula string
        (Composition("Fe2O3"), -10.0, -0.4),  # composition
        ("Fe4O6", -20.0, -0.4),  # complex composition
        ("Fe", -2.0, -1.0),  # single atom
        ("O2", -6.0, -1.0),  # diatomic
    ],
)
def test_calc_energy_from_e_refs_various_inputs(
    input_obj: str | Composition,
    total_energy: float,
    expected: float,
    ref_energies: dict[str, float],
) -> None:
    """Test calculation with various input types and compositions."""
    energy = calc_energy_from_e_refs(input_obj, ref_energies, total_energy)
    assert energy == pytest.approx(expected)


@pytest.mark.parametrize(
    "input_obj,ref_energies,expected",
    [
        (ComputedEntry(Composition("FeO"), -5.0), {"Fe": -1.0, "O": -2.0}, -1.0),
        ({"composition": "FeO", "energy": -5.0}, {"Fe": -1.0, "O": -2.0}, -1.0),
    ],
)
def test_calc_energy_from_e_refs_entry_inputs(
    input_obj: EntryLike,
    ref_energies: dict[str, float],
    expected: float,
) -> None:
    """Test calculation with Entry-like inputs."""
    energy = calc_energy_from_e_refs(input_obj, ref_energies)
    assert energy == pytest.approx(expected)


def test_calc_energy_from_e_refs_error_cases(
    dummy_struct: Structure,
    ref_energies: dict[str, float],
) -> None:
    """Test error handling."""
    # Missing total_energy
    with pytest.raises(ValueError, match="total_energy can't be None"):
        calc_energy_from_e_refs(dummy_struct, ref_energies)

    # Missing reference energy
    with pytest.raises(ValueError, match="Missing reference energies"):
        calc_energy_from_e_refs(dummy_struct, {"Fe": -1.0}, total_energy=-5.0)

    # Invalid input type
    with pytest.raises(TypeError, match="Expected Entry, Structure"):
        calc_energy_from_e_refs([1, 2, 3], ref_energies, total_energy=-5.0)


def test_calc_energy_from_e_refs_equivalence_with_get_e_form(
    ref_energies: dict[str, float],
) -> None:
    """Test equivalence with get_e_form_per_atom."""
    test_entry = ComputedEntry("FeO", -5.0)
    e_form1 = calc_energy_from_e_refs(test_entry, ref_energies)
    e_form2 = get_e_form_per_atom(test_entry, ref_energies)
    assert e_form1 == pytest.approx(e_form2)
