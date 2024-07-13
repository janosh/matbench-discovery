import pandas as pd

from matbench_discovery.data import Key, df_wbm

STRUCT_COL = "sevennet_structure"
E_FORM_COL = "e_form_per_atom_sevennet"

result_csv = "./2024-07-11-sevennet-preds.csv.gz"

csv = pd.read_csv(result_csv)
csv = csv.set_index("material_id")

# NOTE, from code, authors did not used absolute values for this filtering
# It due to the fact that outliers gives very low e_from (both mace & sevennet)
bad_mask = csv[E_FORM_COL] - df_wbm[Key.e_form] < -5
csv[~bad_mask].select_dtypes("number").to_csv(
    "./_2024-07-11-sevennet-preds-no-bad.csv.gz"
)
csv[bad_mask].select_dtypes("number").to_csv("./2024-07-11-sevennet-preds-bad.csv.gz")
