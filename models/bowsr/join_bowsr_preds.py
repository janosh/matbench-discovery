# %%
import os
from glob import glob

import pandas as pd
import pymatviz
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery.data import DataFiles, Model
from matbench_discovery.enums import Task

__author__ = "Janosh Riebesell"
__date__ = "2022-09-22"


# %%
module_dir = os.path.dirname(__file__)
date = "2023-01-20"
energy_model = Model.megnet.label.lower()
glob_pattern = f"{date}-bowsr-{energy_model}-wbm-{Task.IS2RE}/*.json.gz"
file_paths = sorted(glob(f"{module_dir}/{glob_pattern}"))
print(f"Found {len(file_paths):,} files for {glob_pattern = }")

dfs: dict[str, pd.DataFrame] = {}


# %%
for file_path in tqdm(file_paths):
    if file_path in dfs:
        continue
    dfs[file_path] = pd.read_json(file_path).set_index(Key.mat_id)


df_bowsr = pd.concat(dfs.values()).round(4)


# %% compare against WBM formation energy targets to make sure we got sensible results
df_wbm = pd.read_csv(DataFiles.wbm_summary.path).set_index(Key.mat_id)


print(f"{len(df_bowsr) - len(df_wbm)=:,} missing ({len(df_bowsr):,} - {len(df_wbm):,})")


# %% sanity check: since Bowsr uses MEGNet as energy model final BOWSR energy and Megnet
# formation energy should be the same
pymatviz.density_scatter(
    x=df_bowsr.e_form_per_atom_bowsr_megnet,
    y=df_bowsr[f"energy_bowsr_{energy_model}"],
)


# %% remove redundant column after sanity check
df_bowsr = df_bowsr.drop(columns=[f"energy_bowsr_{energy_model}"])


# %%
pymatviz.density_scatter(
    x=df_bowsr.e_form_per_atom_bowsr_megnet,
    y=df_wbm.loc[df_bowsr.index].e_form_per_atom_mp2020_corrected,
)


# %%
out_path = file_paths[0].rsplit("/", 1)[0]
df_bowsr = df_bowsr.round(4)
# save energy and formation energy as fast-loading CSV
df_bowsr.select_dtypes("number").to_csv(f"{out_path}.csv")
df_bowsr.reset_index().to_json(
    f"{out_path}.json.gz", default_handler=lambda x: x.as_dict()
)


# in_path = f"{module_dir}/2023-01-23-bowsr-megnet-wbm-IS2RE.json.gz"
# df_bowsr = pd.read_json(in_path).set_index(Key.mat_id)
