import os

import pandas as pd

from matbench_discovery.data import df_wbm
from matbench_discovery.enums import MbdKey

# uses matbench-discovery matbench-discovery commit ID 012ccfe,
# k_srme commit ID 0269a946, pymatviz v0.15.1
pot_name = "grace_2l_oam_l"


pot_name=pot_name.lower()
e_form_potential_col = f"e_form_per_atom_{pot_name}"
csv_path = f"./results/{pot_name}.csv.gz"
if not os.path.exists(csv_path):
    raise FileNotFoundError(f"Could not find {csv_path=}")

df_preds = pd.read_csv(csv_path).set_index("material_id")

# NOTE this filtering was necessary for both MACE and SevenNet because some outliers
# have extremely low e_form (like -1e40) -- seems not to be necessary for
# well-trained NequIP/Allegro models
bad_mask = abs(df_preds[e_form_potential_col] - df_wbm[MbdKey.e_form_wbm]) > 5
n_preds = len(df_preds[e_form_potential_col].dropna())
print(f"{sum(bad_mask)=} is {sum(bad_mask) / len(df_wbm):.2%} of {n_preds:,}")
df_preds.select_dtypes("number").to_csv(f"./{pot_name}-preds.csv.gz")

filtered_csv_path = f"./{pot_name}-filtered_preds.csv.gz"

df_preds.loc[~bad_mask].select_dtypes("number").to_csv(filtered_csv_path)
