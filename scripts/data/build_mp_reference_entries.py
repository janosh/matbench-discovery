"""Fetch all MP ComputedStructureEntries and derive the MP elemental reference entries
used for formation energies. The DFT ground-truth convex hull is built on the fly from
those entries by matscikit's PhaseDiagram in compile_wbm_test_set.py.
"""

# %%
import gzip
import json
import os

import pandas as pd
import pymatviz as pmv
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import MP_DIR, ROOT, today
from matbench_discovery.data import query_mp_api
from matbench_discovery.energy import (
    Entry,
    formation_energies,
    mp_elemental_ref_energies,
)
from matbench_discovery.enums import DataFiles

os.makedirs(MP_DIR, exist_ok=True)


# %% run on 2022-09-16 and again on 2023-02-07 (then with pymatgen's legacy MPRester)
# one ComputedStructureEntry dict per material: the functional (GGA or GGA+U) the
# GGA/GGA+U mixing scheme picked for it
docs = query_mp_api(
    "materials/thermo", [Key.mat_id, "energy_type", "entries"], thermo_types="GGA_GGA+U"
)

# save all ComputedStructureEntries to disk in the single-document JSON layout that
# DataFiles.mp_computed_structure_entries is read with (columns: material_id, entry)
df_mp_cse = pd.DataFrame(
    {"entry": [doc["entries"][doc["energy_type"]] for doc in docs]},
    index=pd.Index([doc[Key.mat_id] for doc in docs], name=Key.mat_id),
)
# mp-15590 appeared twice on 2023-02-07
df_mp_cse = df_mp_cse[~df_mp_cse.index.duplicated()]
df_mp_cse.reset_index().to_json(
    f"{MP_DIR}/{today}-mp-computed-structure-entries.json.gz"
)


# %%
data_path = DataFiles.mp_computed_structure_entries.path
df_mp_cse = pd.read_json(data_path).set_index(Key.mat_id)

# corrected energies (energy + correction) and compositions, structures not needed
mp_computed_entries = {
    mat_id: Entry.from_dict(dct)
    for mat_id, dct in tqdm(df_mp_cse.entry.items(), total=len(df_mp_cse))
}

print(f"{len(mp_computed_entries)=:,} on {today}")
assert len(mp_computed_entries) == 154_718
# len(mp_computed_entries) = 146,323 on 2022-09-16
# len(mp_computed_entries) = 154,719 on 2023-02-07
# len(mp_computed_entries) = 154,718 on 2024-08-07


# %% compute terminal reference entries across all MP (can be used to compute MP
# compatible formation energies quickly): the lowest-energy entry of each element
elemental_ref_ids: dict[str, str] = {}
for mat_id, entry in mp_computed_entries.items():
    if not entry.composition.is_element:
        continue
    element = entry.composition.elements[0]
    if (
        element not in elemental_ref_ids
        or entry.energy_per_atom
        < mp_computed_entries[elemental_ref_ids[element]].energy_per_atom
    ):
        elemental_ref_ids[element] = mat_id

# save MP phase diagram elemental_ref_entries as JSON
mp_pd_elem_ref_entries_path = f"{MP_DIR}/{today}-mp-elemental-reference-entries.json.gz"
with gzip.open(mp_pd_elem_ref_entries_path, "wt") as file:
    # ComputedEntry layout: the entry dicts without their structures
    json.dump(
        {
            element: {
                key: val
                for key, val in df_mp_cse.entry[mat_id].items()
                if key != "structure"
            }
            | {"@class": "ComputedEntry"}
            for element, mat_id in elemental_ref_ids.items()
        },
        file,
    )


df_mp = pd.read_csv(DataFiles.mp_energies.path, na_filter=False).set_index(Key.mat_id)


# %%
e_form_us = "e_form_us"
mp_entries = [mp_computed_entries[mp_id] for mp_id in df_mp.index]
df_mp[e_form_us] = formation_energies(
    [entry.composition for entry in mp_entries],
    [entry.energy for entry in mp_entries],
    mp_elemental_ref_energies,
)


# make sure get_form_energy_per_atom() reproduces MP formation energies
fig = pmv.density_scatter(df=df_mp, x=Key.formation_energy_per_atom, y=e_form_us)
fig.layout.title = "Formation Energy Comparison"
fig.layout.xaxis.title = "MP Formation Energy (eV/atom)"
fig.layout.yaxis.title = "Our Formation Energy (eV/atom)"
pmv.save_fig(fig, f"{ROOT}/tmp/{today}-our-vs-mp-formation-energies.webp", dpi=300)
