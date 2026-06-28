"""Evaluate MLIP-predicted thermal conductivity metrics against DFT results without
non-analytical correction term (NAC).
"""

# %%
import os

from pymatviz.enums import Key

from matbench_discovery.cli import cli_args
from matbench_discovery.enums import DataFiles
from matbench_discovery.metrics import phonons
from matbench_discovery.phonons import read_kappa_json


def main() -> int:
    """Evaluate kappa metrics and update model YAML files.

    Returns:
        Exit code: 0 if at least one model was evaluated, 1 otherwise.
    """
    models_to_evaluate = cli_args.models
    print(f"Evaluating kappa metrics for {len(models_to_evaluate)} model(s)...")

    n_success = 0
    n_skipped = 0

    for model in models_to_evaluate:
        kappa_103_path = model.kappa_103_path
        if not isinstance(kappa_103_path, str) or not os.path.isfile(kappa_103_path):
            print(f"Skipping {model.label}: no kappa_103_path found")
            n_skipped += 1
            continue

        try:
            print(f"\nProcessing {model.label}...")

            # Load and process data
            df_ml = read_kappa_json(kappa_103_path)
            df_dft = read_kappa_json(DataFiles.phonondb_pbe_103_kappa_no_nac.path)
            df_ml_metrics = phonons.calc_kappa_metrics_from_dfs(df_ml, df_dft)

            # Calculate metrics
            kappa_sre = df_ml_metrics[Key.sre].mean()
            kappa_srme = df_ml_metrics[Key.srme].mean()
            print(f"\t{kappa_srme=:.4f}")
            print(f"\t{kappa_sre=:.4f}")

            # Update YAML file
            metrics_dict = {"srme": kappa_srme, "sre": kappa_sre}
            phonons.write_metrics_to_yaml(model, metrics_dict, kappa_103_path)
            print(f"\tUpdated {model.yaml_path}")
            n_success += 1

        except (ValueError, OSError, KeyError) as exc:
            print(f"\tError processing {model.label}: {exc}")
            n_skipped += 1
            continue

    # Exit with error if no models were successfully evaluated
    if n_success == 0:
        print(f"\nNo models evaluated successfully ({n_skipped} skipped)")
        return 1
    print(f"\nSuccessfully evaluated {n_success} model(s), {n_skipped} skipped")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
