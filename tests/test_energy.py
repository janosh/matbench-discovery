from __future__ import annotations

from typing import Any, Callable

import pytest
from pymatgen.analysis.phase_diagram import PDEntry
from pymatgen.core import Lattice, Structure
from pymatgen.entries.computed_entries import (
    ComputedEntry,
    ComputedStructureEntry,
    Entry,
)

from matbench_discovery.energy import get_e_form_per_atom, get_elemental_ref_entries

dummy_struct = Structure(
    lattice=Lattice.cubic(5),
    species=("Fe", "O"),
    coords=((0, 0, 0), (0.5, 0.5, 0.5)),
)


@pytest.mark.parametrize(
    "constructor", [PDEntry, ComputedEntry, ComputedStructureEntry, lambda **x: x]
)
def test_get_e_form_per_atom(
    constructor: Callable[..., Entry | dict[str, Any]]
) -> None:
    """Test that the formation energy of a composition is computed correctly."""

    entry = {"composition": {"Fe": 1, "O": 1}, "energy": -2.5}
    elemental_ref_entries = {
        "Fe": {"composition": {"Fe": 1}, "energy": -1.0},
        "O": {"composition": {"O": 1}, "energy": -1.0},
    }
    if constructor == ComputedStructureEntry:
        entry["structure"] = dummy_struct
        entry.pop("composition")

    entry = constructor(**entry)

    # don't use ComputedStructureEntry for elemental ref entries, would need many
    # dummy structures
    if constructor == ComputedStructureEntry:
        constructor = ComputedEntry
    elemental_ref_entries = {
        k: constructor(**v) for k, v in elemental_ref_entries.items()
    }
    assert get_e_form_per_atom(entry, elemental_ref_entries) == -0.25


@pytest.mark.parametrize("constructor", [PDEntry, ComputedEntry, lambda **x: x])
@pytest.mark.parametrize("verbose", [True, False])
def test_get_elemental_ref_entries(
    constructor: Callable[..., Entry | dict[str, Any]], verbose: bool
) -> None:
    """Test that the elemental reference entries are correctly identified."""
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
