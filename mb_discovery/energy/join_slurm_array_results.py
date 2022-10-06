# %%
from __future__ import annotations

import os
from datetime import datetime
from glob import glob

import pandas as pd
from tqdm import tqdm

from mb_discovery import ROOT, as_dict_handler

__author__ = "Janosh Riebesell"
__date__ = "2022-08-16"

today = f"{datetime.now():%Y-%m-%d}"


# %%
module_dir = os.path.dirname(__file__)
glob_pattern = f"{module_dir}/2022-10-06-wbm-e-above-hull/*.csv"
file_paths = sorted(glob(glob_pattern))
print(f"Found {len(file_paths):,} files for {glob_pattern = }")

dfs: dict[str, pd.DataFrame] = {}


# %%
for file_path in tqdm(file_paths):
    if file_path in dfs:
        continue
    df = pd.read_csv(file_path).set_index("material_id")
    dfs[file_path] = df


# %%
df_hull = pd.concat(dfs.values())

df_hull.isna().sum()

df_hull_rhys = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wbm-e-above-mp-hull.csv"
).set_index("material_id")

df_hull_rhys.isna().sum()

df_hull["e_above_hull_rhys"] = df_hull_rhys.e_above_mp_hull

ax = df_hull.plot.scatter("e_above_hull_mp_wbm", "e_above_hull_rhys", loglog=False)
ax.axline((0, 0), (1, 1), color="red")
ax.set(xlim=(-1, 1), ylim=(-1, 1))

df_hull_rhys.e_above_mp_hull[df_hull_rhys.e_above_mp_hull.between(-1, 1)].plot.hist(
    bins=100
)
sum(df_hull_rhys.e_above_mp_hull < 0)


# %%
out_path = f"{ROOT}/models/m3gnet/{today}-wbm-e_above_hull.csv"
df_hull.reset_index().to_csv(out_path, default_handler=as_dict_handler)

# out_path = f"{ROOT}/models/m3gnet/2022-08-16-m3gnet-wbm-IS2RE.csv.gz"
# df_m3gnet = pd.read_csv(out_path).set_index("material_id")
