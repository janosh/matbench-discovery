# %%
from __future__ import annotations

import os
from glob import glob

import pandas as pd
import pymatviz
from tqdm import tqdm

from matbench_discovery import ROOT, today

__author__ = "Janosh Riebesell"
__date__ = "2022-09-22"


# %%
module_dir = os.path.dirname(__file__)
task_type = "IS2RE"
date = "2022-11-22"
glob_pattern = f"{date}-bowsr-megnet-wbm-{task_type}/*.json.gz"
file_paths = sorted(glob(f"{module_dir}/{glob_pattern}"))
print(f"Found {len(file_paths):,} files for {glob_pattern = }")

dfs: dict[str, pd.DataFrame] = {}


# %%
for file_path in tqdm(file_paths):
    if file_path in dfs:
        continue
    df = pd.read_json(file_path).set_index("material_id")

    dfs[file_path] = df


# %%
df_bowsr = pd.concat(dfs.values()).round(6)


# %% compare against WBM formation energy targets to make sure we got sensible results
data_path = f"{ROOT}/data/wbm/2022-10-19-wbm-summary.csv"
df_wbm = pd.read_csv(data_path).set_index("material_id")


print(f"{len(df_bowsr):,} - {len(df_wbm):,} = {len(df_bowsr) - len(df_wbm) = :,}")


# %%
pymatviz.density_scatter(
    x=df_bowsr.e_form_per_atom_bowsr_megnet,
    y=df_wbm.loc[df_bowsr.index].e_form_per_atom_mp2020_corrected,
)


# %%
out_path = f"{ROOT}/models/bowsr/{today}-bowsr-megnet-wbm-{task_type}.json.gz"
df_bowsr.reset_index().to_json(out_path, default_handler=lambda x: x.as_dict())

# save energy and formation energy as CSV for fast loading
df_bowsr.select_dtypes("number").to_csv(out_path.replace(".json.gz", ".csv"))

# in_path = f"{ROOT}/models/bowsr/2022-11-22-bowsr-megnet-wbm-IS2RE.json.gz"
# df_bowsr = pd.read_json(in_path).set_index("material_id")
