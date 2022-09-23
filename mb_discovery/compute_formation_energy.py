import itertools

from pymatgen.analysis.phase_diagram import Entry
from tqdm import tqdm


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


def get_form_energy_per_atom(
    entry: Entry, elemental_ref_entries: dict[str, Entry]
) -> float:
    """Get the formation energy of a composition from a list of entries and elemental
    reference energies.
    """
    comp = entry.composition
    form_energy = entry.energy - sum(
        comp[el] * elemental_ref_entries[str(el)].energy_per_atom
        for el in entry.composition.elements
    )

    return form_energy / entry.composition.num_atoms
