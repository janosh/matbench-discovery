"""Formation and elemental reference energies from composition-energy entries."""

import itertools
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import pandas as pd
from matscikit import Composition, Structure
from tqdm import tqdm

from matbench_discovery.enums import DataFiles


@dataclass(frozen=True, init=False)
class Entry:
    """A composition with its total (corrected) energy in eV."""

    composition: Composition
    energy: float
    entry_id: str | None

    def __init__(
        self,
        composition: Composition | str | Mapping[str, float],
        energy: float,
        entry_id: str | None = None,
    ) -> None:
        """Coerce a formula string or element->amount mapping to a Composition."""
        if not isinstance(composition, Composition):
            composition = Composition(composition)
        object.__setattr__(self, "composition", composition)
        object.__setattr__(self, "energy", float(energy))
        object.__setattr__(self, "entry_id", entry_id)

    @property
    def energy_per_atom(self) -> float:
        """Total energy divided by the number of atoms."""
        return self.energy / self.composition.num_atoms

    @classmethod
    def from_dict(cls, dct: Mapping[str, Any]) -> Entry:
        """Build from a serialized entry (pymatgen PDEntry/ComputedEntry layout).

        ComputedEntry dicts store the uncorrected energy plus the total ``correction``
        (the sum of their energy_adjustments), so the corrected energy is
        ``energy + correction`` exactly as pymatgen's ``ComputedEntry.energy``.
        """
        energy = float(dct["energy"]) + float(dct.get("correction", 0.0))
        return cls(dct["composition"], energy, dct.get("entry_id"))


def get_elemental_ref_entries(
    entries: Sequence[Entry | dict[str, Any]], *, verbose: bool = True
) -> dict[str, Entry]:
    """Get the lowest energy Entry for each element in a list of entries.

    Args:
        entries (Sequence[Entry | dict]): Entries or serialized entry dicts to find
            elemental reference entries of.
        verbose (bool, optional): Whether to show a progress bar. Defaults to True.

    Raises:
        ValueError: If some elements are missing terminal reference entries.

    Returns:
        dict[str, Entry]: Map from element symbol to its lowest energy entry.
    """
    normalized_entries: list[Entry] = [
        Entry.from_dict(entry) if isinstance(entry, dict) else entry
        for entry in entries
    ]
    elements = {
        element
        for entry in normalized_entries
        for element in entry.composition.elements
    }
    dim = len(elements)

    if verbose:
        print(
            f"Sorting {len(normalized_entries)} entries with {dim} dimensions...",
            flush=True,
        )

    sorted_entries = sorted(
        normalized_entries, key=lambda entry: entry.composition.reduced_formula
    )

    elemental_ref_entries = {}
    for _, entry_group_iter in tqdm(
        itertools.groupby(
            sorted_entries, key=lambda entry: entry.composition.reduced_formula
        ),
        disable=not verbose,
        desc="Finding elemental reference entries",
    ):
        entry_group = list(entry_group_iter)
        min_entry = min(entry_group, key=lambda entry: entry.energy_per_atom)
        composition = min_entry.composition
        if composition.is_element:
            elemental_ref_entries[composition.elements[0]] = min_entry

    element_symbols = set(map(str, elements))
    if len(elemental_ref_entries) < dim:
        missing = element_symbols - set(elemental_ref_entries)
        raise ValueError(f"Missing terminal entries: {missing}")

    return elemental_ref_entries


# contains all MP elemental reference entries to compute formation energies
# produced by get_elemental_ref_entries() in build_phase_diagram.py
mp_elem_ref_entries = (
    pd.read_json(DataFiles.mp_elemental_ref_entries.path, typ="series")
    .map(Entry.from_dict)
    .to_dict()
)

# tested to agree with TRI's MP reference energies
# https://github.com/TRI-AMDD/CAMD/blob/1c965cba636/camd/utils/data.py#L134
mp_elemental_ref_energies = {
    elem: entry.energy_per_atom for elem, entry in mp_elem_ref_entries.items()
}


def calc_energy_from_e_refs(
    struct_or_entry: Entry | dict[str, Any] | Structure | Composition | str,
    ref_energies: dict[str, float],
    total_energy: float | None = None,
) -> float:
    """Calculate energy per atom relative to reference states (e.g., for formation or
    cohesive energy calculations).

    Args:
        struct_or_entry (Entry | Structure | Composition | str): Either:
            - An Entry or entry dict containing 'energy' and 'composition' keys
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

        if isinstance(struct_or_entry, str):
            comp = Composition(struct_or_entry)
        elif isinstance(struct_or_entry, Structure):
            comp = struct_or_entry.composition
        elif isinstance(struct_or_entry, Composition):
            comp = struct_or_entry
        else:
            cls_name = type(struct_or_entry).__name__
            raise TypeError(
                "Expected Entry, Structure, Composition or formula string, "
                f"got {cls_name}"
            )

    if missing_refs := set(comp.elements) - set(ref_energies):
        raise ValueError(f"Missing reference energies for elements: {missing_refs}")

    e_ref = sum(ref_energies[el] * amt for el, amt in comp.items())

    return (energy - e_ref) / comp.num_atoms
