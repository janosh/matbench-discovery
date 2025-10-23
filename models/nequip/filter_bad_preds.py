import os

import pandas as pd

from matbench_discovery.data import df_wbm
from matbench_discovery.enums import MbdKey

# uses commits matbench-discovery 012ccfe, k_srme commit 0269a946, pymatviz v0.15.1

e_form_nequip_col = "e_form_per_atom_nequip"
csv_path = "./results/nequip.csv.gz"
if not os.path.isfile(csv_path):
    csv_path = "./nequip.csv.gz"

df_preds = pd.read_csv(csv_path).set_index("material_id")

# NOTE this filtering was necessary for both MACE and SevenNet because some outliers
# have extremely low e_form (like -1e40) -- seems not to be necessary for
# well-trained NequIP/Allegro models
bad_mask = abs(df_preds[e_form_nequip_col] - df_wbm[MbdKey.e_form_wbm]) > 5
n_preds = len(df_preds[e_form_nequip_col].dropna())
print(f"{sum(bad_mask)=} is {sum(bad_mask) / len(df_wbm):.2%} of {n_preds:,}")
df_preds.select_dtypes("number").to_csv("./nequip-preds.csv.gz")
df_preds.loc[~bad_mask].select_dtypes("number").to_csv("./nequip-filtered_preds.csv.gz")
