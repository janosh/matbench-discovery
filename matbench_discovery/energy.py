import itertools
from collections.abc import Sequence

import pandas as pd
from pymatgen.analysis.phase_diagram import Entry, PDEntry
from pymatgen.core import Composition, Structure
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry
from pymatgen.util.typing import EntryLike
from tqdm import tqdm

from matbench_discovery.data import DATA_FILES


def get_elemental_ref_entries(
    entries: Sequence[EntryLike], verbose: bool = True
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
mp_elem_reference_entries = (
    pd.read_json(DATA_FILES.mp_elemental_ref_entries, typ="series")
    .map(ComputedEntry.from_dict)
    .to_dict()
)

# tested to agree with TRI's MP reference energies
# https://github.com/TRI-AMDD/CAMD/blob/1c965cba636531e542f4821a555b98b2d81ed034/camd/utils/data.py#L134
mp_elemental_ref_energies = {
    elem: round(entry.energy_per_atom, 4)
    for elem, entry in mp_elem_reference_entries.items()
}


def get_e_form_per_atom(
    entry: EntryLike,
    elemental_ref_energies: dict[str, float] = mp_elemental_ref_energies,
) -> float:
    """Get the formation energy of a composition from a list of entries and a dict
    mapping elements to reference energies.

    Args:
        entry: Entry | dict[str, float | str | Composition]: pymatgen Entry (PDEntry,
            ComputedEntry or ComputedStructureEntry) or dict with energy and composition
            keys to compute formation energy of.
        elemental_ref_energies (dict[str, float], optional): Must be a covering set (for
            entry) of terminal reference energies, i.e. eV/atom of the lowest energy
            elemental phase for each element. Defaults to MP elemental reference
            energies as collected on 2022-09-19 get_elemental_ref_entries(). This was
            tested to give the same formation energies as found in MP.

    Returns:
        float: formation energy in eV/atom.
    """
    if isinstance(entry, dict):
        energy = entry["energy"]
        comp = Composition(entry["composition"])  # is idempotent if already Composition
    elif isinstance(entry, Entry):
        energy = entry.energy
        comp = entry.composition
    else:
        raise TypeError(
            f"{entry=} must be Entry (or subclass like ComputedEntry) or dict"
        )

    e_refs = {str(el): elemental_ref_energies[str(el)] for el in comp}

    for key, ref_entry in e_refs.items():
        if isinstance(ref_entry, dict):
            e_refs[key] = PDEntry.from_dict(ref_entry)

    form_energy = energy - sum(comp[el] * e_refs[str(el)] for el in comp)

    return form_energy / comp.num_atoms


df_cse: pd.DataFrame = None


def apply_mp_2020_corrections_to_ml_struct(
    df_ml: pd.DataFrame,
    energy_col: str,
    struct_col: str,
    entry_col: str = "computed_structure_entry",
) -> pd.DataFrame:
    """MaterialsProject2020Compatibility are structure-dependent for oxides and
    sulfides. This function applies MP2020 corrections to a DataFrame of
    ComputedStructureEntries where the DFT structure is replaced by one predicted by a
    ML force field.

    Args:
        df_ml (pd.DataFrame): DataFrame with ML-predicted energies and structures.
        energy_col (str): Column name of ML-predicted energies.
        struct_col (str): Column name of ML-predicted structures.
        entry_col (str, optional): Column name of ComputedStructureEntries. Defaults to
            "computed_structure_entry".

    Returns:
        pd.DataFrame: DataFrame with MP2020-corrected energies and structures.
    """
    id_col = "material_id"

    # assign to global df_cse so data is only loaded once per python session
    global df_cse  # noqa: PLW0603
    if df_cse is None:
        print(
            "Loading WBM ComputedStructureEntries. This takes a bit on the 1st call "
            "to this function",
            flush=True,
        )
        df_cse = pd.read_json(DATA_FILES.wbm_computed_structure_entries)
        df_cse = df_cse.set_index(id_col)
        df_cse[entry_col] = [
            ComputedStructureEntry.from_dict(dct) for dct in tqdm(df_cse[entry_col])
        ]

    # assign ML-predicted energies and ML-relaxed structures to WBM
    # ComputedStructureEntries. needed because MP2020 energy
    # corrections applied below are structure-dependent (for oxides and sulfides)
    cse: ComputedStructureEntry
    for row in tqdm(df_ml.itertuples(), total=len(df_ml)):
        mat_id = row.Index
        ml_energy = row[energy_col]
        mlip_struct = Structure.from_dict(row[struct_col])
        df_ml.loc[mat_id, struct_col] = mlip_struct
        cse = df_cse.loc[mat_id, entry_col]
        cse._energy = ml_energy  # cse._energy is the uncorrected energy
        cse._structure = mlip_struct
        df_ml.loc[mat_id, entry_col] = cse

    # apply energy corrections
    processed = MaterialsProject2020Compatibility().process_entries(
        df_ml[entry_col], verbose=True, clean=True
    )
    assert len(processed) == len(df_ml)  # make sure no entries were skipped

    return df_ml
