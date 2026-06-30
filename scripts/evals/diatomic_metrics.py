"""Calculate diatomic curve metrics for all models and write them to YAML files."""

import gzip
import json
import numbers
import os
import sys

from matbench_discovery import ROOT
from matbench_discovery.cli import cli_args
from matbench_discovery.enums import MbdKey
from matbench_discovery.metrics import diatomics
from matbench_discovery.metrics.diatomics import DiatomicCurves
from matbench_discovery.remote.fetch import maybe_auto_download_file
from models.run_diatomics import DIATOMIC_METRIC_EXCLUSIONS, drop_metric_exclusions

models_to_evaluate = cli_args.models
print(f"Evaluating diatomic metrics for {len(models_to_evaluate)} model(s)...")

n_success = 0
n_skipped = 0

metrics_to_write: dict[str, dict[str, object]] = {
    metric: {}
    for metric in (
        MbdKey.tortuosity,
        MbdKey.energy_diff_flips,
        MbdKey.energy_jump,
        MbdKey.force_flips,
        MbdKey.force_total_variation,
        MbdKey.force_jump,
        MbdKey.pbe_wall_dist_mae,
        MbdKey.pbe_energy_mae,
        MbdKey.pbe_bond_length_error,
        MbdKey.pbe_well_depth_error,
        MbdKey.pbe_force_mae,
        MbdKey.pbe_vib_freq_error,
    )
}

pbe_ref_curves = diatomics.load_dft_reference_curves("PBE")

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

    metrics = diatomics.calc_diatomic_metrics(
        ref_curves=pbe_ref_curves,
        pred_curves=pred_curves,
        metrics=metrics_to_write,
        interpolate=200,
    )
    metrics = drop_metric_exclusions(model.name, metrics)
    excluded_formulas = list(DIATOMIC_METRIC_EXCLUSIONS.get(model.name, ()))

    # Write metrics to YAML
    mean_metrics = diatomics.write_metrics_to_yaml(
        model,
        metrics,
        run_metadata={"excluded_formulas": excluded_formulas},
    )
    print(f"{model.label}:")
    for metric, val in mean_metrics.items():
        value_str = f"{val:.5}" if isinstance(val, numbers.Real) else str(val)
        print(f"  {metric}: {value_str}")
    n_success += 1

# Exit with error if no models were successfully evaluated
if n_success == 0:
    print(f"\nNo models evaluated successfully ({n_skipped} skipped)")
    sys.exit(1)
print(f"\nSuccessfully evaluated {n_success} model(s), {n_skipped} skipped")
