"""Calculate diatomic curve metrics for all models and write them to YAML files."""

import gzip
import json
import os
import sys

from matbench_discovery import ROOT
from matbench_discovery.cli import cli_args
from matbench_discovery.enums import Model
from matbench_discovery.metrics import diatomics
from matbench_discovery.metrics.diatomics import DiatomicCurves
from matbench_discovery.remote.fetch import maybe_auto_download_file

models_to_evaluate = cli_args.models or list(Model)
print(f"Evaluating diatomic metrics for {len(models_to_evaluate)} model(s)...")

n_success = 0
n_skipped = 0

# Loop over all models
for model in models_to_evaluate:
    if not os.path.isfile(model.yaml_path):
        print(f"Skipping {model.label}: YAML file not found")
        n_skipped += 1
        continue

    diatomics_metrics = model.metrics.get("diatomics")
    if not isinstance(diatomics_metrics, dict):
        print(f"Skipping {model.label}: no diatomics metrics config in YAML")
        n_skipped += 1
        continue

    pred_file = diatomics_metrics.get("pred_file")
    if not isinstance(pred_file, str):
        print(f"Skipping {model.label}: no pred_file in diatomics config")
        n_skipped += 1
        continue

    abs_path = f"{ROOT}/{pred_file}"
    pred_file_url = diatomics_metrics.get("pred_file_url")

    # Try to download if file doesn't exist and URL is available
    if not os.path.isfile(abs_path) and pred_file_url:
        maybe_auto_download_file(
            pred_file_url, abs_path, label=f"{model.label} diatomics"
        )

    if not os.path.isfile(abs_path):
        print(f"Skipping {model.label}: prediction file not found at {pred_file}")
        n_skipped += 1
        continue

    # Load predicted curves
    with gzip.open(abs_path, mode="rb") as file:
        pred_data = json.load(file) or {}

    pred_curves = DiatomicCurves.from_dict(pred_data)

    # Calculate metrics (without reference data)
    metrics = diatomics.calc_diatomic_metrics(ref_curves=None, pred_curves=pred_curves)

    # Write metrics to YAML
    mean_metrics = diatomics.write_metrics_to_yaml(model, metrics)
    print(f"✓ {model.label}:")
    for metric, val in mean_metrics.items():
        print(f"  {metric}: {val:.5}")
    n_success += 1

# Exit with error if no models were successfully evaluated
if n_success == 0:
    print(f"\n✗ No models evaluated successfully ({n_skipped} skipped)")
    sys.exit(1)
else:
    print(f"\n✓ Successfully evaluated {n_success} model(s), {n_skipped} skipped")
