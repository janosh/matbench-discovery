"""
Templated from nequip/join_matbench_preds.py nequip/extract_wbm_structs.py
"""

# uses commits matbench-discovery f0e54b7

from glob import glob

import pandas as pd
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery.data import as_dict_handler

e_form_tace_col = "e_form_per_atom_tace"
struct_col = "tace_structure"
results = "results"
pot_name = "tace"
out_path = f"{results}/{pot_name}"
files = sorted(glob(f"{results}/{pot_name}-*.json.gz"))

dfs = {}
for file_path in tqdm(files, desc="Loading results"):
    if file_path in dfs:
        continue
    df_i = pd.read_json(file_path).set_index(Key.mat_id)
    dfs[file_path] = df_i

df_tace = pd.concat(dfs.values()).round(4)
df_tace.reset_index().to_json(
    "tace-wbm-geo-opt.jsonl.gz",
    lines=True,
    default_handler=as_dict_handler,
    orient="records",
)
