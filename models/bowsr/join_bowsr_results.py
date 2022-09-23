# %%
from __future__ import annotations

import os
from datetime import datetime
from glob import glob

import pandas as pd
from pymatgen.core import Structure
from tqdm import tqdm

from mb_discovery import ROOT, as_dict_handler
from mb_discovery.plots import hist_classified_stable_as_func_of_hull_dist

__author__ = "Janosh Riebesell"
__date__ = "2022-09-22"

today = f"{datetime.now():%Y-%m-%d}"


# %%
module_dir = os.path.dirname(__file__)
task_type = "IS2RE"
date = "2022-09-22"
glob_pattern = f"{date}-bowsr-wbm-{task_type}/*.json.gz"
file_paths = sorted(glob(f"{module_dir}/{glob_pattern}"))
print(f"Found {len(file_paths):,} files for {glob_pattern = }")

dfs: dict[str, pd.DataFrame] = {}


# %%
# 2022-08-16 tried multiprocessing.Pool() to load files in parallel but was somehow
# slower than serial loading
for file_path in tqdm(file_paths):
    if file_path in dfs:
        continue
    # keep whole dataframe in memory
    df = pd.read_json(file_path).set_index("material_id")
    col_map = dict(
        structure_pred="structure_bowsr",
        energy_pred="energy_bowsr",
        e_form_per_atom_pred="e_form_per_atom_bowsr",
    )
    df = df.rename(columns=col_map)
    df["structure_bowsr"] = df.structure_bowsr.map(Structure.from_dict)
    df["formula"] = df.structure_bowsr.map(lambda x: x.formula)
    df["volume"] = df.structure_bowsr.map(lambda x: x.volume)
    df["n_sites"] = df.structure_bowsr.map(len)
    dfs[file_path] = df


# %%
df_bowsr = pd.concat(dfs.values())


# %%
df_wbm = pd.read_csv(  # download wbm-steps-summary.csv (23.31 MB)
    "https://figshare.com/files/37570234?private_link=ff0ad14505f9624f0c05"
).set_index("material_id")

df_bowsr["e_form_wbm"] = df_wbm.e_form_per_atom


# %%
df_bowsr.hist(bins=200, figsize=(18, 12))
df_bowsr.isna().sum()


# %%
out_path = f"{ROOT}/models/bowsr/{today}-bowsr-wbm-{task_type}.json.gz"
df_bowsr.reset_index().to_json(out_path, default_handler=as_dict_handler)

out_path = f"{ROOT}/models/bowsr/2022-08-16-bowsr-wbm-IS2RE.json.gz"
df_bowsr = pd.read_json(out_path).set_index("material_id")


# %%
df_hull = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wbm-e-above-mp-hull.csv"
).set_index("material_id")
df_bowsr["e_above_mp_hull"] = df_hull.e_above_mp_hull
df_bowsr["e_above_hull_pred"] = (  # TODO fix this incorrect e_above_hull_pred
    df_bowsr["e_form_per_atom_bowsr"] - df_bowsr["e_above_mp_hull"]
)

ax_hull_dist_hist = hist_classified_stable_as_func_of_hull_dist(
    e_above_hull_pred=df_bowsr.e_above_hull_pred,
    e_above_hull_true=df_bowsr.e_above_mp_hull,
)

# ax_hull_dist_hist.figure.savefig(f"{ROOT}/plots/{today}-bowsr-wbm-hull-dist-hist.pdf")
