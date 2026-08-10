"""Calculate diatomic curve metrics for all models and write them to YAML files."""

import gzip
import json
import numbers
import os

from matbench_discovery import ROOT
from matbench_discovery.cli import cli_args
from matbench_discovery.data import file_ref_name, file_ref_url
from matbench_discovery.enums import MbdKey, Model
from matbench_discovery.metrics import diatomics
from matbench_discovery.metrics.diatomics import DiatomicCurves
from matbench_discovery.metrics.diatomics.exclusions import drop_metric_exclusions
from matbench_discovery.remote.fetch import maybe_auto_download_file
from scripts.evals import evaluate_models

METRICS_TO_WRITE: dict[str, dict[str, object]] = {
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


def main() -> int:
    """Evaluate diatomic metrics, returning 0 if any model succeeds and 1 otherwise."""
    pbe_ref_curves = diatomics.load_dft_reference_curves("PBE")

    def evaluate_one(model: Model) -> str | None:
        if not os.path.isfile(model.yaml_path):
            return "YAML file not found"

        diatomics_metrics = model.metrics.get("diatomics") or {}
        if not diatomics_metrics:
            return "no diatomics metrics config in YAML"

        pred_file = diatomics_metrics.get("pred_file")
        pred_name = file_ref_name(pred_file)
        if pred_name is None:
            return "no pred_file in diatomics config"

        prediction_path = f"{ROOT}/{pred_name}"
        pred_file_url = file_ref_url(pred_file)

        if not os.path.isfile(prediction_path) and pred_file_url:
            maybe_auto_download_file(
                pred_file_url, prediction_path, label=f"{model.label} diatomics"
            )

        if not os.path.isfile(prediction_path):
            return f"prediction file not found at {pred_name}"

        with gzip.open(prediction_path, mode="rb") as file:
            pred_data = json.load(file) or {}

        metrics = diatomics.calc_diatomic_metrics(
            ref_curves=pbe_ref_curves,
            pred_curves=DiatomicCurves.from_dict(pred_data),
            metrics=METRICS_TO_WRITE,
            interpolate=200,
        )
        metrics = drop_metric_exclusions(model.name, metrics)
        # Preserve source metadata by leaving run_metadata unset.
        mean_metrics = diatomics.write_metrics_to_yaml(model, metrics)
        print(f"{model.label}:")
        for metric, value in mean_metrics.items():
            value_str = f"{value:.5}" if isinstance(value, numbers.Real) else str(value)
            print(f"  {metric}: {value_str}")
        return None

    return evaluate_models("diatomic", cli_args.models, evaluate_one)


if __name__ == "__main__":
    raise SystemExit(main())
