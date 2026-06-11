"""Aggregate MD metric prediction files and write them to model YAML files."""

from __future__ import annotations

import os
import sys

import pandas as pd

from matbench_discovery.cli import cli_args
from matbench_discovery.enums import Model
from matbench_discovery.metrics import md


def main() -> int:
    """Evaluate available MD prediction files for selected models."""
    models_to_evaluate = cli_args.models or list(Model)
    print(f"Evaluating MD metrics for {len(models_to_evaluate)} model(s)...")

    n_success = 0
    n_skipped = 0

    for model in models_to_evaluate:
        md_path = model.md_path
        if not md_path or not os.path.isfile(md_path):
            print(f"Skipping {model.label}: no metrics.md.pred_file found")
            n_skipped += 1
            continue

        try:
            print(f"\nProcessing {model.label}...")
            df_eval = pd.read_csv(md_path)
            metrics = md.calc_md_metrics(df_eval)
            md_config = model.metrics.get("md")
            pred_file_url = (
                md_config.get("pred_file_url") if isinstance(md_config, dict) else None
            )
            md.write_metrics_to_yaml(
                model,
                metrics,
                md_path,
                pred_file_url=(
                    pred_file_url if isinstance(pred_file_url, str) else None
                ),
            )
            for key, value in metrics.items():
                print(f"\t{key}={value:.4f}")
            print(f"\tUpdated {model.yaml_path}")
            n_success += 1
        except Exception as exc:  # noqa: BLE001 - continue evaluating other models
            print(f"\tError processing {model.label}: {exc}")
            n_skipped += 1

    if n_success == 0:
        print(f"\nNo models evaluated successfully ({n_skipped} skipped)")
        return 1
    print(f"\nSuccessfully evaluated {n_success} model(s), {n_skipped} skipped")
    return 0


if __name__ == "__main__":
    sys.exit(main())
