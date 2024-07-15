import pandas as pd

from matbench_discovery.data import Key, df_wbm

E_FORM_COL = "e_form_per_atom_sevennet"

csv_path = "./2024-07-11-sevennet-preds.csv.gz"
df_preds = pd.read_csv(csv_path).set_index("material_id")

# NOTE this filtering is necessary for both MACE and SevenNet because some outliers
# have extremely low e_form (like -1e40)
bad_mask = df_preds[E_FORM_COL] - df_wbm[Key.e_form] < -5
df_preds[~bad_mask].select_dtypes("number").to_csv(
    "./2024-07-11-sevennet-preds-no-bad.csv.gz"
)
df_preds[bad_mask].select_dtypes("number").to_csv(
    "./2024-07-11-sevennet-preds-bad.csv.gz"
)
