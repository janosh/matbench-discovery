import pandas as pd

from matbench_discovery import ROOT

df_wbm = pd.read_csv(f"{ROOT}/data/wbm/2022-10-19-wbm-summary.csv")
df_wbm = df_wbm.set_index("material_id")
