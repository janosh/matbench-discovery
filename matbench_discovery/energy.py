"""Functions to compute formation and elemental reference energies from
pymatgen EntryLikes.
"""

from collections.abc import Sequence
from typing import Any

import pandas as pd
from pymatgen.analysis.phase_diagram import Entry, PDEntry
from pymatgen.core import Composition, Structure
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.util.typing import EntryLike
from tqdm import tqdm

from matbench_discovery.enums import DataFiles


def get_elemental_ref_entries(
    entries: Sequence[EntryLike | dict[str, Any]], *, verbose: bool = True
) -> dict[str, Entry]:
    """Get the lowest energy pymatgen Entry for each element in a list of entries.

    Args:
        entries (Sequence[Entry]): pymatgen Entries (PDEntry, ComputedEntry or
            ComputedStructureEntry) to find elemental reference entries of.
        verbose (bool, optional): Whether to show a progress bar. Defaults to True.

    Raises:
        ValueError: If some elements are missing terminal reference entries.

    Returns:
        dict[str, Entry]: Map from element symbol to its lowest energy entry.
    """
    elemental_ref_entries: dict[str, Entry] = {}
    elements: set[str] = set()
    for entry_like in tqdm(
        entries, disable=not verbose, desc="Finding elemental reference entries"
    ):
        entry = (
            PDEntry.from_dict(entry_like)
            if isinstance(entry_like, dict)
            else entry_like
        )
        comp = entry.composition
        elements.update(map(str, comp))
        if comp.is_element:
            symbol = str(comp.elements[0])
            ref_entry = elemental_ref_entries.get(symbol)
            if ref_entry is None or entry.energy_per_atom < ref_entry.energy_per_atom:
                elemental_ref_entries[symbol] = entry

    if missing := elements - set(elemental_ref_entries):
        raise ValueError(f"Missing terminal entries: {missing}")

    return elemental_ref_entries


# contains all MP elemental reference entries to compute formation energies
# produced by get_elemental_ref_entries() in build_mp_reference_entries.py
mp_elem_ref_entries = (
    pd.read_json(DataFiles.mp_elemental_ref_entries.path, typ="series")
    .map(ComputedEntry.from_dict)
    .to_dict()
)

# tested to agree with TRI's MP reference energies
# https://github.com/TRI-AMDD/CAMD/blob/1c965cba636/camd/utils/data.py#L134
mp_elemental_ref_energies = {
    elem: entry.energy_per_atom for elem, entry in mp_elem_ref_entries.items()
}


def calc_energy_from_e_refs(
    struct_or_entry: EntryLike | dict[str, Any] | Structure | Composition | str,
    ref_energies: dict[str, float],
    total_energy: float | None = None,
) -> float:
    """Calculate energy per atom relative to reference states (e.g., for formation or
    cohesive energy calculations).

    Args:
        struct_or_entry (EntryLike | Structure | Composition | str): Either:
            - A pymatgen Entry (PDEntry, ComputedEntry, etc.) or entry dict containing
              'energy' and 'composition' keys
            - A Structure or Composition object or formula string (must also provide
              total_energy)
        ref_energies (dict[str, float]): Dictionary of reference energies per atom.
            For formation energy: elemental reference energies (e.g.
            mp_elemental_ref_energies).
            For cohesive energy: isolated atom reference energies
        total_energy (float | None): Total energy of the structure/composition. Required
            if struct_or_entry is not an Entry or entry dict. Ignored otherwise.

    Returns:
        float: Energy per atom relative to references (e.g., formation or cohesive
        energy) in the same units as input energies.

    Raises:
        TypeError: If input types are invalid
        ValueError: If missing reference energies for some elements
    """
    if isinstance(struct_or_entry, dict):
        if missing_keys := {"composition", "energy"} - struct_or_entry.keys():
            raise ValueError(
                f"Entry dict missing required keys: {sorted(missing_keys)}"
            )
        energy: float = struct_or_entry["energy"]
        comp = Composition(struct_or_entry["composition"])
    elif isinstance(struct_or_entry, Entry):
        energy = struct_or_entry.energy
        comp = struct_or_entry.composition
    else:
        if total_energy is None:
            raise ValueError("total_energy can't be None when 1st arg is not an Entry")
        energy = total_energy

        if isinstance(struct_or_entry, Structure):
            comp = struct_or_entry.composition
        elif isinstance(struct_or_entry, str | Composition):
            comp = Composition(struct_or_entry)
        else:
            cls_name = type(struct_or_entry).__name__
            raise TypeError(
                "Expected Entry, Structure, Composition or formula string, "
                f"got {cls_name}"
            )

    if missing_refs := set(map(str, comp)) - set(ref_energies):
        raise ValueError(f"Missing reference energies for elements: {missing_refs}")

    e_ref = sum(ref_energies[str(el)] * amt for el, amt in comp.items())

    return (energy - e_ref) / comp.num_atoms
