"""Calculate diatomic curve metrics for all models and write them to YAML files."""

import gzip
import json
import numbers

from matbench_discovery.cli import cli_args
from matbench_discovery.enums import Model
from matbench_discovery.metrics import diatomics
from matbench_discovery.metrics.diatomics import DiatomicCurves
from scripts.evals import evaluate_models


def main() -> int:
    """Evaluate diatomic metrics with the shared sweep's failure and skip handling."""
    pbe_ref_curves = diatomics.load_dft_reference_curves("PBE")

    def evaluate_one(model: Model) -> str | None:
        """Score one model's diatomic predictions and write its metrics to YAML."""
        # getter validates the artifact path and downloads/md5-verifies the file
        if (prediction_path := model.diatomics_path) is None:
            return "no pred_file in diatomics config"

        with gzip.open(prediction_path, mode="rb") as file:
            pred_data = json.load(file)

        pred_curves = DiatomicCurves.from_dict(pred_data)
        metrics = diatomics.score_diatomic_curves(
            model.name, pbe_ref_curves, pred_curves
        )
        mean_metrics = diatomics.write_metrics_to_yaml(
            model,
            metrics,
            # YAML cost is re-derived from the file's own run_metadata
            run_metadata=pred_data.get("run_metadata", {}),
            vib_freq_coverage=diatomics.calc_vib_freq_coverage(
                pbe_ref_curves, pred_curves, metrics
            ),
        )
        print(f"{model.label}:")
        for metric, value in mean_metrics.items():
            value_str = f"{value:.5}" if isinstance(value, numbers.Real) else str(value)
            print(f"  {metric}: {value_str}")
        return None

    return evaluate_models("diatomic", cli_args.models, evaluate_one)


if __name__ == "__main__":
    raise SystemExit(main())
