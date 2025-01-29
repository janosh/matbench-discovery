"""Evaluate MLIP-predicted thermal conductivity metrics against DFT results without
non-analytical correction term (NAC).
"""

# %%
import pandas as pd
from pymatviz.enums import Key

from matbench_discovery.data import DataFiles, Model
from matbench_discovery.metrics import phonons

model_name = "mace-omat-0-medium"
in_folder = f"2025-01-28-{model_name}-phononDB-LTC-FIRE_2SR_force0.0001_sym1e-05"
model_dir = Model.mace_mpa_0.yaml_path.rsplit("/", 1)[0]
in_path = f"{model_dir}/{in_folder}/kappa.json.gz"

df_ml = pd.read_json(in_path).set_index(Key.mat_id)

json_path = DataFiles.phonondb_pbe_103_kappa_no_nac.path
df_dft = pd.read_json(json_path).set_index("mp_id")

df_ml_metrics = phonons.calc_kappa_metrics_from_dfs(df_ml, df_dft)

kappa_sre = df_ml_metrics[Key.sre].mean()
kappa_srme = df_ml_metrics[Key.srme].mean()

print(f"{model_name = }\n\t{kappa_srme = :.6f}\n\t{kappa_sre = :.6f}")
