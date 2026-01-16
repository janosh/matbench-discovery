"""Evaluate MLIP-predicted thermal conductivity metrics against DFT results without
non-analytical correction term (NAC).
"""

# %%
import os
import sys

import pandas as pd
from pymatviz.enums import Key

from matbench_discovery.cli import cli_args
from matbench_discovery.enums import DataFiles, Model
from matbench_discovery.metrics import phonons


def main() -> int:
    """Evaluate kappa metrics and update model YAML files.

    Returns:
        Exit code: 0 if at least one model was evaluated, 1 otherwise.
    """
    models_to_evaluate = cli_args.models or list(Model)
    print(f"Evaluating kappa metrics for {len(models_to_evaluate)} model(s)...")

    n_success = 0
    n_skipped = 0

    for model in models_to_evaluate:
        if not os.path.isfile(model.kappa_103_path or ""):
            print(f"Skipping {model.label}: no kappa_103_path found")
            n_skipped += 1
            continue

        try:
            print(f"\nProcessing {model.label}...")

            # Load and process data
            df_ml = pd.read_json(model.kappa_103_path).set_index(Key.mat_id)
            df_dft = pd.read_json(DataFiles.phonondb_pbe_103_kappa_no_nac.path)
            if "mp_id" in df_dft.columns:
                df_dft = df_dft.rename(columns={"mp_id": Key.mat_id})
            df_dft = df_dft.set_index(Key.mat_id)
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
            n_success += 1

        except Exception as exc:
            print(f"\t✗ Error processing {model.label}: {exc}")
            n_skipped += 1
            continue

    # Exit with error if no models were successfully evaluated
    if n_success == 0:
        print(f"\n✗ No models evaluated successfully ({n_skipped} skipped)")
        return 1
    print(f"\n✓ Successfully evaluated {n_success} model(s), {n_skipped} skipped")
    return 0


if __name__ == "__main__":
    sys.exit(main())
