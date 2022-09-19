# %%
import json
from datetime import datetime

import pandas as pd
from aviary.wren.utils import get_aflow_label_from_spglib
from pymatgen.analysis.phase_diagram import PDEntry
from pymatgen.entries.computed_entries import ComputedStructureEntry
from tqdm import tqdm

from mb_discovery import ROOT, as_dict_handler
from mb_discovery.compute_formation_energy import get_form_energy_per_atom

"""
Change JSON orientation of wbm-cleaned.json.gz and WBM step IDs to match the dielectric
Pareto frontier project.
"""

__author__ = "Janosh Riebesell"
__date__ = "2022-06-26"

today = f"{datetime.now():%Y-%m-%d}"


def increment_wbm_material_id(wbm_id: str) -> str:
    """Turns 'step_1-0' into 'wbm-step-1-1' etc."""
    *_, step_num, material_num = wbm_id.split("_")

    assert step_num.isdigit() and material_num.isdigit()

    return f"wbm-step-{step_num}-{int(material_num) + 1}"


# %%
df_wbm = pd.read_json(
    f"{ROOT}/data/2022-06-11-from-rhys/wbm-cleaned.json.gz", orient="split"
).set_index("material_id")
df_wbm.index = df_wbm.index.map(increment_wbm_material_id)

df_wbm = df_wbm.rename(columns=dict(bandgap="bandgap_pbe"))

df_wbm.cse = [ComputedStructureEntry.from_dict(x) for x in tqdm(df_wbm.cse)]

for wbm_id, cse in tqdm(df_wbm.cse.items()):  # type: ignore
    cse.entry_id = wbm_id

df_wbm["energy"] = [cse.energy for cse in df_wbm.cse]
df_wbm["n_sites"] = [len(cse.structure) for cse in df_wbm.cse]
df_wbm["formula"] = [cse.structure.formula for cse in df_wbm.cse]


df_wbm.reset_index().to_json(
    f"{ROOT}/data/{today}-wbm-cses-and-initial-structures.json.gz",
    default_handler=as_dict_handler,
)


# %% compute WBM formation energies and Aflow-style Wyckoff labels
# warning: slow (takes ~15 min)
df_wbm["wyckoff"] = [
    get_aflow_label_from_spglib(cse.structure) for cse in tqdm(df_wbm.cse)
]

df_wren = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wren-mp-initial-structures.csv"
).set_index("material_id")
df_wren.e_form_target.head(10)
# load MP elemental reference entries
with open(f"{ROOT}/data/2022-09-19-mp-elemental-reference-entries.json") as json_file:
    mp_elemental_reference_entries = json.load(json_file)
    for elem, entry in mp_elemental_reference_entries.items():
        mp_elemental_reference_entries[elem] = PDEntry.from_dict(entry)

df_wbm["e_form_per_atom_wrt_mp_elemental_refs"] = [
    get_form_energy_per_atom(entry, mp_elemental_reference_entries)
    for entry in tqdm(df_wbm.cse)
]

df_wbm.drop(columns=["initial_structure", "cse"]).to_csv(
    f"{ROOT}/data/{today}-wbm-formation-energy+wyckoff.csv"
)

max_e_form = 3
# the file 2022-09-20-wbm-formation-energy+wyckoff was uploaded to figshare as
# wbm-steps-summary.csv (23.5 MB)
# https://figshare.com/files/37542841
df_wbm.drop(columns=["initial_structure", "cse"], errors="ignore").query(
    f"e_form_per_atom < {max_e_form}"
).to_csv(f"{ROOT}/data/{today}-wbm-formation-energy+wyckoff.csv")

# df_wbm = pd.read_csv(
#     f"{ROOT}/data/2022-09-19-wbm-formation-energy+wyckoff.csv"
# ).set_index("material_id")


# %% 2022-07-18 also increment material_ids in wbm-e-above-mp-hull.csv
for filename in tqdm(
    "wbm-e-above-mp-hull",
    "wren-mp-initial-structures",
    "cgcnn-mp-initial-structures",
    "voronoi-mp-initial-structures",
    "wren-mp-cse",
    "cgcnn-mp-cse",
    "voronoi-mp-cse",
):
    file_path = f"{ROOT}/data/2022-06-11-from-rhys/{filename}.csv"

    df = pd.read_csv(file_path).round(5)

    df["material_id"] = df.material_id.map(increment_wbm_material_id)

    df.to_csv(file_path, index=False)
