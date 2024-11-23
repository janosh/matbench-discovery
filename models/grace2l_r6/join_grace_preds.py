import glob
import os

import pandas as pd
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery.data import as_dict_handler, df_wbm
from matbench_discovery.energy import calc_energy_from_e_refs, mp_elemental_ref_energies
from matbench_discovery.enums import MbdKey, Task

__author__ = "Yury Lysogorskiy"
__date__ = "2024-11-22"

energy_column = "grace2l_r6_energy"
e_form_grace_col = "e_form_per_atom_grace"
struct_col = "grace2l_r6_structure"


module_dir = os.path.dirname(__file__)
task_type = Task.IS2RE
date = "2024-11-21"
glob_pattern = "2024-11-21-MP_GRACE_2L_r6_11Nov2024-wbm-IS2RE-FIRE/production-*.json.gz"
file_paths = glob.glob(glob_pattern)

print(f"Found {len(file_paths):,} files for {glob_pattern = }")


dfs: list[pd.DataFrame] = []
for fn in file_paths:
    print(fn)
    dfs.append(pd.read_json(fn))


tot_df = pd.concat(dfs)
tot_df["id_tuple"] = (
    tot_df["material_id"].str.split("-").map(lambda x: (int(x[1]), int(x[2])))
)
tot_df = (
    tot_df.sort_values("id_tuple")
    .reset_index(drop=True)
    .drop(columns=["id_tuple", struct_col])
)

df_grace = tot_df.set_index("material_id")
df_grace[Key.formula] = df_wbm[Key.formula]


print("Calculating formation energies")
e_form_list = []
for _, row in tqdm(df_grace.iterrows(), total=len(df_grace)):
    e_form = calc_energy_from_e_refs(
        row["formula"],
        ref_energies=mp_elemental_ref_energies,
        total_energy=row[energy_column],
    )
    e_form_list.append(e_form)


df_grace[e_form_grace_col] = e_form_list

df_wbm[[*df_grace]] = df_grace


# %%
bad_mask = abs(df_wbm[e_form_grace_col] - df_wbm[MbdKey.e_form_dft]) > 5
n_preds = len(df_wbm[e_form_grace_col].dropna())
print(f"{sum(bad_mask)=} is {sum(bad_mask) / len(df_wbm):.2%} of {n_preds:,}")
out_path = file_paths[0].rsplit("/", 1)[0]

df_grace = df_grace.round(4)
df_grace.select_dtypes("number").to_csv(f"{out_path}.csv.gz")


df_grace.reset_index().to_json(f"{out_path}.json.gz", default_handler=as_dict_handler)

df_bad = df_grace[bad_mask].copy()
df_bad[MbdKey.e_form_dft] = df_wbm[MbdKey.e_form_dft]
df_bad.to_csv(f"{out_path}-bad.csv")
