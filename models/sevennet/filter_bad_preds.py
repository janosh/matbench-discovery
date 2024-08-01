import pandas as pd

from matbench_discovery.data import df_wbm
from matbench_discovery.enums import MbdKey

e_form_7net_col = "e_form_per_atom_sevennet"

csv_path = "./2024-07-11-sevennet-preds.csv.gz"
df_preds = pd.read_csv(csv_path).set_index("material_id")

# NOTE this filtering is necessary for both MACE and SevenNet because some outliers
# have extremely low e_form (like -1e40)
bad_mask = abs(df_preds[e_form_7net_col] - df_wbm[MbdKey.e_form_wbm]) > 5
n_preds = len(df_preds[e_form_7net_col].dropna())
print(f"{sum(bad_mask)=} is {sum(bad_mask) / len(df_wbm):.2%} of {n_preds:,}")
df_preds.select_dtypes("number").to_csv("./2024-07-11-sevennet-preds.csv.gz")
