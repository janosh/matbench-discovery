"""Evaluate MLIP-predicted thermal conductivity metrics against DFT results without
non-analytical correction term (NAC).
"""

# %%
import os

import pandas as pd
from pymatviz.enums import Key

from matbench_discovery.data import DataFiles, Model
from matbench_discovery.metrics import phonons

for model in Model:
    if model.kappa_103_path is None or not os.path.isfile(model.kappa_103_path):
        continue

    df_ml = pd.read_json(model.kappa_103_path).set_index(Key.mat_id)

    json_path = DataFiles.phonondb_pbe_103_kappa_no_nac.path
    df_dft = pd.read_json(json_path).set_index(Key.mat_id)

    df_ml_metrics = phonons.calc_kappa_metrics_from_dfs(df_ml, df_dft)

    kappa_sre = df_ml_metrics[Key.sre].mean()
    kappa_srme = df_ml_metrics[Key.srme].mean()

    print(f"{model.label=}\n\t{kappa_srme=:.4f}\n\t{kappa_sre=:.4f}")
