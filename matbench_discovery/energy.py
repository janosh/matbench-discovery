"""Functions to compute formation and elemental reference energies from
pymatgen EntryLikes.
"""

import itertools
import warnings
from collections.abc import Sequence
from typing import Any

import pandas as pd
from pymatgen.analysis.phase_diagram import Entry, PDEntry
from pymatgen.core import Composition, Structure
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.util.typing import EntryLike
from tqdm import tqdm

from matbench_discovery.data import DataFiles


def get_elemental_ref_entries(
    entries: Sequence[EntryLike], *, verbose: bool = True
) -> dict[str, Entry]:
    """Get the lowest energy pymatgen Entry for each element in a list of entries.

    Args:
        entries (Sequence[Entry]): pymatgen Entries (PDEntry, ComputedEntry or
            ComputedStructureEntry) to find elemental reference entries of.
        verbose (bool, optional): Whether to show a progress bar. Defaults to False.

    Raises:
        ValueError: If some elements are missing terminal reference entries.
        ValueError: If there are more terminal entries than dimensions. Should never
            happen.

    Returns:
        dict[str, Entry]: Map from element symbol to its lowest energy entry.
    """
    entries = [PDEntry.from_dict(e) if isinstance(e, dict) else e for e in entries]
    elements = {elems for entry in entries for elems in entry.composition.elements}
    dim = len(elements)

    if verbose:
        print(f"Sorting {len(entries)} entries with {dim} dimensions...", flush=True)

    entries = sorted(entries, key=lambda e: e.composition.reduced_composition)

    elemental_ref_entries = {}
    for composition, entry_group in tqdm(
        itertools.groupby(entries, key=lambda e: e.composition.reduced_composition),
        disable=not verbose,
        desc="Finding elemental reference entries",
    ):
        min_entry = min(entry_group, key=lambda e: e.energy_per_atom)
        if composition.is_element:
            elem_symb = str(composition.elements[0])
            elemental_ref_entries[elem_symb] = min_entry

    if len(elemental_ref_entries) > dim:
        missing = elements - set(elemental_ref_entries)
        raise ValueError(f"Some terminal entries are {missing = }")
    if len(elemental_ref_entries) < dim:
        extra = set(elemental_ref_entries) - set(elements)
        raise ValueError(
            f"There are more terminal element entries than dimensions: {extra}"
        )

    return elemental_ref_entries


# contains all MP elemental reference entries to compute formation energies
# produced by get_elemental_ref_entries() in build_phase_diagram.py
mp_elem_ref_entries = (
    pd.read_json(DataFiles.mp_elemental_ref_entries.path, typ="series")
    .map(ComputedEntry.from_dict)
    .to_dict()
)

# tested to agree with TRI's MP reference energies
# https://github.com/TRI-AMDD/CAMD/blob/1c965cba636531e542f4821a555b98b2d81ed034/camd/utils/data.py#L134
mp_elemental_ref_energies = {
    elem: round(entry.energy_per_atom, 4) for elem, entry in mp_elem_ref_entries.items()
}


def calc_energy_from_e_refs(
    struct_or_entry: EntryLike | Structure | Composition | str,
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
    if isinstance(struct_or_entry, dict):  # entry dict case
        energy = struct_or_entry["energy"]
        comp = Composition(struct_or_entry["composition"])
    # Entry/ComputedEntry/ComputedStructureEntry instance case
    elif isinstance(struct_or_entry, Entry):
        energy = struct_or_entry.energy
        comp = struct_or_entry.composition
    else:  # Structure/Composition/formula case
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

    # Check that we have all needed reference energies
    if missing_refs := set(map(str, comp)) - set(ref_energies):
        raise ValueError(f"Missing reference energies for elements: {missing_refs}")

    # Calculate reference energy
    e_ref = sum(ref_energies[str(el)] * amt for el, amt in comp.items())

    return (energy - e_ref) / comp.num_atoms


def get_e_form_per_atom(*args: Any, **kwargs: Any) -> float:  # noqa: D417
    """Get formation energy for a phase diagram entry (1st arg, composition + absolute
    energy) and a dict mapping elements to per-atom reference energies (2nd arg).

    Args:
        entry: Entry | dict[str, float | str | Composition]: pymatgen Entry (PDEntry,
            ComputedEntry or ComputedStructureEntry) or dict with energy (absolute, not
            per atom) and composition keys to compute formation energy of.
        elemental_ref_energies (dict[str, float], optional): Must be a covering set (for
            entry) of terminal reference energies, i.e. eV/atom of the lowest energy
            elemental phase for each element. Defaults to MP elemental reference
            energies as collected on 2022-09-19 get_elemental_ref_entries(). This was
            tested to give the same formation energies as found in MP.

    Returns:
        float: formation energy in eV/atom.

    Raises:
        TypeError: If entry is not a pymatgen Entry or dict.
    """
    warnings.warn(
        "get_e_form_per_atom is deprecated. Use calc_energy_from_e_refs instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    if len(args) >= 2:
        ref_energies = args[1]
        args = (args[0], *args[2:])
    else:
        if entry := kwargs.pop("entry"):
            args = (entry, *args)
        ref_energies = kwargs.pop("elemental_ref_energies", mp_elemental_ref_energies)
    kwargs.setdefault("ref_energies", ref_energies)
    return calc_energy_from_e_refs(*args, **kwargs)
