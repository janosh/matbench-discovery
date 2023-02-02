# %%
from __future__ import annotations

import os
import warnings
from glob import glob

import pandas as pd
from pymatgen.core import Structure
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatviz import density_scatter
from tqdm import tqdm

from matbench_discovery import ROOT, today
from matbench_discovery.data import as_dict_handler
from matbench_discovery.energy import get_e_form_per_atom

__author__ = "Janosh Riebesell"
__date__ = "2022-08-16"

warnings.filterwarnings(action="ignore", category=UserWarning, module="pymatgen")


# %%
module_dir = os.path.dirname(__file__)
task_type = "IS2RE"
date = "2022-10-31"
glob_pattern = f"{date}-m3gnet-wbm-{task_type}/*.json.gz"
file_paths = sorted(glob(f"{module_dir}/{glob_pattern}"))
print(f"Found {len(file_paths):,} files for {glob_pattern = }")

dfs: dict[str, pd.DataFrame] = {}


# %%
for file_path in tqdm(file_paths):
    if file_path in dfs:
        continue
    df = pd.read_json(file_path).set_index("material_id")
    df[f"m3gnet_energy_{task_type}"] = [
        x["energies"][-1][0] for x in df.m3gnet_trajectory
    ]
    # drop trajectory to save memory
    dfs[file_path] = df.drop(columns="m3gnet_trajectory")


# %%
df_m3gnet = pd.concat(dfs.values()).round(4)


# %%
df_wbm = pd.read_json(
    f"{ROOT}/data/wbm/2022-10-19-wbm-computed-structure-entries.json.bz2"
).set_index("material_id")

df_wbm["cse"] = [
    ComputedStructureEntry.from_dict(x) for x in tqdm(df_wbm.computed_structure_entry)
]


# %% transfer M3GNet energies and relaxed structures WBM CSEs
# make sure we're not skipping m3gnet structures, other way around is fine
# i.e. WBM CSEs for which M3GNet structures are missing
assert not {*df_m3gnet.index} - {*df_wbm.index}

cse: ComputedStructureEntry
for cse in tqdm(df_wbm.cse):
    if cse.entry_id in df_m3gnet.index:
        # cse._energy is the uncorrected energy
        cse._energy = df_m3gnet.loc[cse.entry_id, "m3gnet_energy"]
        m3gnet_struct = Structure.from_dict(
            df_m3gnet.loc[cse.entry_id, "m3gnet_structure"]
        )
        cse._structure = m3gnet_struct


# %%
df_wbm["e_form_per_atom_m3gnet_uncorrected"] = [
    get_e_form_per_atom(cse) for cse in tqdm(df_wbm.cse)
]


# %% apply energy corrections
out = MaterialsProject2020Compatibility().process_entries(
    df_wbm.cse, verbose=True, clean=True
)
assert len(out) == len(df_wbm)


# %% compute corrected formation energies
df_wbm["e_form_per_atom_m3gnet"] = [
    get_e_form_per_atom(cse) for cse in tqdm(df_wbm.cse)
]

e_form_m3gnet_cols = list(df_wbm.filter(like="m3gnet"))
assert len(e_form_m3gnet_cols) == 2
df_m3gnet[e_form_m3gnet_cols] = df_wbm[e_form_m3gnet_cols]


# %%
ax = density_scatter(
    df=df_m3gnet, x="e_form_per_atom_m3gnet", y="e_form_per_atom_m3gnet_uncorrected"
)


# %%
out_path = f"{module_dir}/{today}-m3gnet-wbm-{task_type}.json.gz"
df_m3gnet = df_m3gnet.round(4)
df_m3gnet.reset_index().to_json(out_path, default_handler=as_dict_handler)

df_m3gnet.select_dtypes("number").to_csv(out_path.replace(".json.gz", ".csv"))

# in_path = f"{module_dir}/2022-10-31-m3gnet-wbm-IS2RE.json.gz"
# df_m3gnet = pd.read_json(in_path).set_index("material_id")
