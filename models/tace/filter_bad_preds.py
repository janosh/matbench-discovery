"""
Templated from nequip/filter_bad_preds.py
"""

# uses commits matbench-discovery f0e54b7

import os

import pandas as pd

from matbench_discovery.data import df_wbm
from matbench_discovery.enums import MbdKey

e_form_tace_col = "e_form_per_atom_tace"
csv_path = "results/tace.csv.gz"
if not os.path.isfile(csv_path):
    csv_path = "tace.csv.gz"

df_preds = pd.read_csv(csv_path).set_index("material_id")

# NOTE this filtering was necessary for both MACE and SevenNet because some outliers
# have extremely low e_form (like -1e40) -- seems not to be necessary for
# well-trained tace/Allegro models
# TACE follow them to filter
bad_mask = abs(df_preds[e_form_tace_col] - df_wbm[MbdKey.e_form_wbm]) > 5
n_preds = len(df_preds[e_form_tace_col].dropna())
print(f"{sum(bad_mask)=} is {sum(bad_mask) / len(df_wbm):.2%} of {n_preds:,}")
df_preds.select_dtypes("number").to_csv("tace-preds.csv.gz")
df_preds.loc[~bad_mask].select_dtypes("number").to_csv("tace-filtered_preds.csv.gz")
