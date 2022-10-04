import itertools

import pandas as pd
from pymatgen.analysis.phase_diagram import Entry
from pymatgen.entries.computed_entries import ComputedEntry
from tqdm import tqdm

from mb_discovery import ROOT


def get_elemental_ref_entries(
    entries: list[Entry], verbose: bool = False
) -> dict[str, Entry]:

    elements = {elems for entry in entries for elems in entry.composition.elements}
    dim = len(elements)

    if verbose:
        print(f"Sorting {len(entries)} entries with {dim} dimensions...")
    entries = sorted(entries, key=lambda e: e.composition.reduced_composition)

    elemental_ref_entries = {}
    if verbose:
        print("Finding elemental reference entries...", flush=True)
    for composition, group in tqdm(
        itertools.groupby(entries, key=lambda e: e.composition.reduced_composition)
    ):
        min_entry = min(group, key=lambda e: e.energy_per_atom)
        if composition.is_element:
            elem_symb = str(composition.elements[0])
            elemental_ref_entries[elem_symb] = min_entry

    if len(elemental_ref_entries) > dim:
        missing = elements - set(elemental_ref_entries)
        raise ValueError(f"Some terminal entries are {missing = }")
    elif len(elemental_ref_entries) < dim:
        extra = set(elemental_ref_entries) - set(elements)
        raise ValueError(
            f"There are more terminal element entries than dimensions: {extra}"
        )

    return elemental_ref_entries


# contains all MP elemental reference entries to compute formation energies
# produced by get_elemental_ref_entries() in build_phase_diagram.py
mp_elem_refs_path = f"{ROOT}/data/2022-09-19-mp-elemental-reference-entries.json"
try:
    mp_elem_reference_entries = (
        pd.read_json(mp_elem_refs_path, typ="series")
        .map(ComputedEntry.from_dict)
        .to_dict()
    )
except FileNotFoundError:
    mp_elem_reference_entries = None


def get_form_energy_per_atom(
    entry: Entry, elemental_ref_entries: dict[str, Entry] = None
) -> float:
    """Get the formation energy of a composition from a list of entries and elemental
    reference energies.

    Args:
        entry (Entry): pymatgen Entry (PDEntry, ComputedEntry or ComputedStructureEntry)
            to compute formation energy of.
        elemental_ref_entries (dict[str, Entry], optional): Must be a complete set of
            terminal (i.e. elemental) reference entries containing the lowest energy
            phase for each element present in entry. Defaults to MP elemental reference
            entries as collected on 2022-09-19 get_elemental_ref_entries(). This was
            tested to give the same formation energies as computed by MP.

    Returns:
        float: formation energy in eV/atom.
    """
    if elemental_ref_entries is None:
        if mp_elem_reference_entries is None:
            raise ValueError(
                f"Couldn't load {mp_elem_refs_path=}, you must pass "
                f"{elemental_ref_entries=} explicitly."
            )

        elemental_ref_entries = mp_elem_reference_entries

    comp = entry.composition
    form_energy = entry.energy - sum(
        comp[el] * elemental_ref_entries[str(el)].energy_per_atom
        for el in entry.composition.elements
    )

    return form_energy / entry.composition.num_atoms
