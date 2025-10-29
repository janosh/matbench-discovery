"""Evaluate MLIP-predicted thermal conductivity metrics against DFT results without
non-analytical correction term (NAC).
"""

# %%
import os

import pandas as pd
from pymatviz.enums import Key

from matbench_discovery.cli import cli_args
from matbench_discovery.enums import DataFiles, Model
from matbench_discovery.metrics import phonons


def main() -> None:
    """Evaluate kappa metrics and update model YAML files."""
    models_to_evaluate = [Model.matris_10m_oam, Model.matris_10m_mp]
    for model in models_to_evaluate:

        if not os.path.isfile(model.kappa_103_path or ""):
            print(f"Skipping {model.label}: no kappa_103_path found")
            continue
        
        try:
            print(f"\nProcessing {model.label}...")

            # Load and process data
            df_ml = pd.read_json(model.kappa_103_path).set_index(Key.mat_id)
            """
                * TODO: Bug report
                * version: pymatviz == 0.17.2
                * Key.mat_id = "material_id" while "mp_id" in json file
            """
            df_dft = pd.read_json(
                DataFiles.phonondb_pbe_103_kappa_no_nac.path
            ).set_index("mp_id") # Key.mat_id 
            df_ml_metrics = phonons.calc_kappa_metrics_from_dfs(df_ml, df_dft)

            # Calculate metrics
            kappa_sre = df_ml_metrics[Key.sre].mean()
            kappa_srme = df_ml_metrics[Key.srme].mean()
            print(f"\t{kappa_srme=:.4f}")
            print(f"\t{kappa_sre=:.4f}")

            # Update YAML file
            metrics_dict = {"srme": kappa_srme, "sre": kappa_sre}
            phonons.write_metrics_to_yaml(model, metrics_dict, model.kappa_103_path)
            print(f"\t✓ Updated {model.yaml_path}")

        except Exception as exc:
            print(f"\t✗ Error processing {model.label}: {exc}")
            continue


if __name__ == "__main__":
    main()
