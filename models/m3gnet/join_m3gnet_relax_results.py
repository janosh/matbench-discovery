# %%
from __future__ import annotations

import os
from datetime import datetime
from glob import glob

import pandas as pd
from pymatgen.analysis.phase_diagram import PDEntry
from pymatgen.core import Structure
from tqdm import tqdm

from mb_discovery import ROOT, as_dict_handler
from mb_discovery.energy import get_form_energy_per_atom

__author__ = "Janosh Riebesell"
__date__ = "2022-08-16"

today = f"{datetime.now():%Y-%m-%d}"


# %%
module_dir = os.path.dirname(__file__)
task_type = "RS2RE"
date = "2022-08-19"
glob_pattern = f"{date}-m3gnet-wbm-relax-{task_type}/*.json.gz"
file_paths = sorted(glob(f"{module_dir}/{glob_pattern}"))
print(f"Found {len(file_paths):,} files for {glob_pattern = }")

dfs: dict[str, pd.DataFrame] = {}


# %%
# 2022-08-16 tried multiprocessing.Pool() to load files in parallel but was somehow
# slower than serial loading
for file_path in tqdm(file_paths):
    if file_path in dfs:
        continue
    try:
        # keep whole dataframe in memory
        df = pd.read_json(file_path)
        df.index = df.index.str.replace("_", "-")
        df.index.name = "material_id"
        col_map = dict(
            final_structure="m3gnet_structure", trajectory="m3gnet_trajectory"
        )
        df = df.rename(columns=col_map)
        df.reset_index().to_json(file_path)
        df["m3gnet_energy"] = df.m3gnet_trajectory.map(lambda x: x["energies"][-1][0])
        df["m3gnet_structure"] = df.m3gnet_structure.map(Structure.from_dict)
        df["formula"] = df.m3gnet_structure.map(lambda x: x.formula)
        df["volume"] = df.m3gnet_structure.map(lambda x: x.volume)
        df["n_sites"] = df.m3gnet_structure.map(len)
        dfs[file_path] = df.drop(columns=["m3gnet_trajectory"])
    except (ValueError, FileNotFoundError):
        # pandas v1.5+ correctly raises FileNotFoundError, below raises ValueError
        continue


# %%
df_m3gnet = pd.concat(dfs.values())
if any(df_m3gnet.index.str.contains("_")):
    df_m3gnet.index = df_m3gnet.index.str.replace("_", "-")


# %%
pd_entries_m3gnet = [
    PDEntry(row.m3gnet_structure.composition, row.m3gnet_energy)
    for row in df_m3gnet.itertuples()
]
df_m3gnet["e_form_m3gnet_from_ppd"] = [
    get_form_energy_per_atom(entry) for entry in pd_entries_m3gnet
]


# %% compare against WBM formation energy targets to make sure we got sensible results
df_hull = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wbm-e-above-mp-hull.csv"
).set_index("material_id")

df_m3gnet["e_above_hull_mp"] = df_hull.e_above_hull_mp


# %%
df_m3gnet.hist(bins=200, figsize=(18, 12))
df_m3gnet.isna().sum()


# %%
out_path = f"{ROOT}/models/m3gnet/{today}-m3gnet-wbm-relax-{task_type}.json.gz"
df_m3gnet.reset_index().to_json(out_path, default_handler=as_dict_handler)

# out_path = f"{ROOT}/models/m3gnet/2022-08-16-m3gnet-wbm-IS2RE.json.gz"
# df_m3gnet = pd.read_json(out_path).set_index("material_id")
