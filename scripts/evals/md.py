"""Aggregate per-system MD metric files and write model-level metrics to YAML."""

import os
import sys

import pandas as pd

from matbench_discovery.cli import cli_args
from matbench_discovery.metrics import md as md_metrics


def main() -> int:
    """Evaluate MD metrics and update model YAML files.

    Returns:
        Exit code: 0 if at least one model was evaluated, 1 otherwise.
    """
    models_to_evaluate = cli_args.models
    print(f"Evaluating MD metrics for {len(models_to_evaluate)} model(s)...")

    n_success = 0
    n_skipped = 0

    for model in models_to_evaluate:
        md_path = model.md_path
        if not isinstance(md_path, str) or not os.path.isfile(md_path):
            print(f"Skipping {model.label}: no md_path found")
            n_skipped += 1
            continue

        try:
            print(f"\nProcessing {model.label}...")
            df_md = pd.read_csv(md_path)
            metrics = md_metrics.calc_md_metrics(df_md)
            for key, value in metrics.items():
                print(f"\t{key}={value:.4f}")

            md_yaml = model.metrics.get("md") or {}
            md_metrics.write_metrics_to_yaml(
                model,
                metrics,
                pred_file_path=md_yaml.get("pred_file") or md_path,
                pred_file_url=md_yaml.get("pred_file_url"),
            )
            print(f"\tUpdated {model.yaml_path}")
            n_success += 1
        except (ValueError, OSError, KeyError) as exc:
            print(f"\tError processing {model.label}: {exc}")
            n_skipped += 1

    if n_success == 0:
        print(f"\nNo models evaluated successfully ({n_skipped} skipped)")
        return 1
    print(f"\nSuccessfully evaluated {n_success} model(s), {n_skipped} skipped")
    return 0


if __name__ == "__main__":
    sys.exit(main())
