"""
This processing script has been copied from the 7net script here:
https://github.com/janosh/matbench-discovery/blob/main/models/sevennet/join_7net_preds.py
And then slightly refactored for NequIP/Allegro, and changing the WBM missing structures
error to a warning.
Note that it requires pymatviz >=0.15.0

Takes about 4.5 mins to run.
"""

# uses commits matbench-discovery 012ccfe, k_srme commit 0269a946, pymatviz v0.15.1

from glob import glob

import pandas as pd
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery.data import as_dict_handler

e_form_allegro_col = "e_form_per_atom_allegro"
struct_col = "allegro_structure"
results = "./results"
pot_name = "allegro"
out_path = f"{results}/{pot_name}"
files = sorted(glob(f"{results}/{pot_name}-*.json.gz"))

dfs = {}
for file_path in tqdm(files, desc="Loading results"):
    if file_path in dfs:
        continue
    df_i = pd.read_json(file_path).set_index(Key.mat_id)
    dfs[file_path] = df_i

df_allegro = pd.concat(dfs.values()).round(4)
df_allegro.reset_index().to_json(
    "allegro-wbm-geo-opt.jsonl.gz",
    lines=True,
    default_handler=as_dict_handler,
    orient="records",
)
