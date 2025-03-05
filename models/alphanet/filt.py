import pandas as pd

from matbench_discovery.data import df_wbm
from matbench_discovery.enums import MbdKey

e_form_anet_col = "e_form_per_atom_alphanet"

csv_path = "./res_relax/alphanet.csv.gz"
df_preds = pd.read_csv(csv_path).set_index("material_id")


bad_mask = abs(df_preds[e_form_anet_col] - df_wbm[MbdKey.e_form_wbm]) > 5
n_preds = len(df_preds[e_form_anet_col].dropna())
print(f"{sum(bad_mask)=} is {sum(bad_mask) / len(df_wbm):.2%} of {n_preds:,}")
df_preds.select_dtypes("number").to_csv("./res_relax/alphanet-no-bad.csv.gz")
