"""Fetch all MP ComputedStructureEntries and derive the MP elemental reference entries
used for formation energies. The DFT ground-truth convex hull is built on the fly from
those entries by a pymatgen PatchedPhaseDiagram in compile_wbm_test_set.py.
"""

# %%
import gzip
import json
import os

import pandas as pd
import pymatviz as pmv
from pymatgen.entries.computed_entries import ComputedEntry
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import MP_DIR, ROOT, today
from matbench_discovery.data import query_mp_api
from matbench_discovery.energy import (
    calc_energy_from_e_refs,
    get_elemental_ref_entries,
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

# drop the structure, just load ComputedEntry (corrected energies and compositions)
mp_computed_entries = [ComputedEntry.from_dict(dct) for dct in tqdm(df_mp_cse.entry)]

print(f"{len(mp_computed_entries)=:,} on {today}")
assert len(mp_computed_entries) == 154_718
# len(mp_computed_entries) = 146,323 on 2022-09-16
# len(mp_computed_entries) = 154,719 on 2023-02-07
# len(mp_computed_entries) = 154,718 on 2024-08-07


# %% compute terminal reference entries across all MP (can be used to compute MP
# compatible formation energies quickly)
elemental_ref_entries = get_elemental_ref_entries(mp_computed_entries)

# save MP phase diagram elemental_ref_entries as JSON
mp_pd_elem_ref_entries_path = f"{MP_DIR}/{today}-mp-elemental-reference-entries.json.gz"
with gzip.open(mp_pd_elem_ref_entries_path, "wt") as file:
    json.dump(elemental_ref_entries, file, default=lambda x: x.as_dict())


df_mp = pd.read_csv(DataFiles.mp_energies.path, na_filter=False).set_index(Key.mat_id)


# %%
e_form_us = "e_form_us"
# entry IDs carry a functional suffix (mp-111-GGA), so key by the material ID
mp_entries_by_id = dict(zip(df_mp_cse.index, mp_computed_entries, strict=True))
df_mp[e_form_us] = [
    calc_energy_from_e_refs(mp_entries_by_id[mp_id], mp_elemental_ref_energies)
    for mp_id in df_mp.index
]


# make sure get_form_energy_per_atom() reproduces MP formation energies
fig = pmv.density_scatter(df=df_mp, x=Key.formation_energy_per_atom, y=e_form_us)
fig.layout.title = "Formation Energy Comparison"
fig.layout.xaxis.title = "MP Formation Energy (eV/atom)"
fig.layout.yaxis.title = "Our Formation Energy (eV/atom)"
pmv.save_fig(fig, f"{ROOT}/tmp/{today}-our-vs-mp-formation-energies.webp", dpi=300)
