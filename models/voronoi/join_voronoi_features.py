# %%
from __future__ import annotations

import os
from glob import glob

import pandas as pd
from tqdm import tqdm

__author__ = "Janosh Riebesell"
__date__ = "2022-08-16"


# %%
module_dir = os.path.dirname(__file__)
date = "2022-11-18"
glob_pattern = f"{date}-voronoi-features-wbm/voronoi-features-wbm-*.csv.bz2"
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
        df = pd.read_csv(file_path).set_index("material_id")
        dfs[file_path] = df
    except FileNotFoundError:
        print(f"{file_path=} not found")
        continue


# %%
df_features = pd.concat(dfs.values())

assert df_features.isna().sum().max() <= 18


# %%
out_path = f"{module_dir}/{date}-voronoi-features-wbm.csv.bz2"
df_features.to_csv(out_path)
