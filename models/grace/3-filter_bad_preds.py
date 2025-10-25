"""Filter outlier predictions (>5 eV/atom deviation from WBM reference)."""

import os

import pandas as pd

from matbench_discovery.data import df_wbm
from matbench_discovery.enums import MbdKey

model_name = "grace_2l_oam_l"
csv_path = f"./results/{model_name}.csv.gz"

if not os.path.isfile(csv_path):
    raise FileNotFoundError(f"{csv_path=}")

df_preds = pd.read_csv(csv_path).set_index("material_id")

bad_mask = (
    abs(df_preds[f"e_form_per_atom_{model_name}"] - df_wbm[MbdKey.e_form_wbm]) > 5
)
n_bad = sum(bad_mask)
n_total = len(df_wbm)
print(f"Filtered {n_bad} ({n_bad / n_total:.2%}) outlier predictions")

df_preds.select_dtypes("number").to_csv(f"./{model_name}-preds.csv.gz")
df_preds.loc[~bad_mask].select_dtypes("number").to_csv(
    f"./{model_name}-filtered_preds.csv.gz"
)
