# %%
import gzip
import itertools
import json
import os
import pickle
from datetime import datetime

import pandas as pd
from pymatgen.analysis.phase_diagram import PatchedPhaseDiagram, PDEntry
from pymatgen.ext.matproj import MPRester
from tqdm import tqdm

from mb_discovery import ROOT

today = f"{datetime.now():%Y-%m-%d}"
module_dir = os.path.dirname(__file__)


# %%
def get_elemental_ref_entries(
    entries: list[PDEntry], verbose: bool = False
) -> dict[str, PDEntry]:

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
    entry: PDEntry, elemental_ref_entries: dict[str, PDEntry]
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


# %%
if __name__ == "__main__":
    all_mp_entries = MPRester().get_entries("")  # run on 2022-09-16
    # mp-15590 appears twice so we drop_duplicates()
    df_mp_entries = pd.DataFrame(all_mp_entries, columns=["entry"]).drop_duplicates()
    df_mp_entries["material_id"] = [x.entry_id for x in df_mp_entries.entry]
    df_mp_entries = df_mp_entries.set_index("material_id")

    df_mp_entries.reset_index().to_json(
        f"{ROOT}/data/{today}-2-all-mp-entries.json.gz",
        default_handler=lambda x: x.as_dict(),
    )

    df_mp_entries = pd.read_json(
        f"{ROOT}/data/2022-09-16-all-mp-entries.json.gz"
    ).set_index("material_id")
    all_mp_entries = [PDEntry.from_dict(x) for x in df_mp_entries.entry]

    print(f"{len(df_mp_entries) = :,}")
    # len(df_mp_entries) = 146,323

    ppd_mp = PatchedPhaseDiagram(all_mp_entries)
    # prints:
    # PatchedPhaseDiagram
    #   Covering 44805 Sub-Spaces

    # save MP PPD to disk
    with gzip.open(f"{module_dir}/{today}-ppd-mp.pkl.gz", "wb") as zip_file:
        pickle.dump(ppd_mp, zip_file)

    elemental_ref_entries = get_elemental_ref_entries(all_mp_entries)

    # save elemental_ref_entries to disk as json
    with open(f"{module_dir}/{today}-elemental-ref-entries.json", "w") as f:
        json.dump(elemental_ref_entries, f, default=lambda x: x.as_dict())
