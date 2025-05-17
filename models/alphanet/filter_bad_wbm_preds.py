"""Filter out bad predictions from AlphaNet on the WBM test set."""

import os

import pandas as pd

from matbench_discovery.data import df_wbm
from matbench_discovery.enums import MbdKey

module_dir = os.path.dirname(__file__)
e_form_anet_col = "e_form_per_atom_alphanet"

csv_path = f"{module_dir}/res_relax/alphanet.csv.gz"
df_preds = pd.read_csv(csv_path).set_index("material_id")

# formation energies that differ by more than 5 eV/atom
bad_mask = abs(df_preds[e_form_anet_col] - df_wbm[MbdKey.e_form_wbm]) > 5
n_preds = len(df_preds[e_form_anet_col].dropna())
print(f"{sum(bad_mask)=} is {sum(bad_mask) / len(df_wbm):.2%} of {n_preds:,}")
csv_no_bad_path = f"{module_dir}/res_relax/alphanet-no-bad.csv.gz"
df_preds.select_dtypes("number").to_csv(csv_no_bad_path)
