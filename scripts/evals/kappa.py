"""Evaluate MLIP-predicted thermal conductivity metrics against DFT results without
non-analytical correction term (NAC).
"""

import os

from matbench_discovery.cli import cli_args
from matbench_discovery.enums import Model
from matbench_discovery.metrics import phonons
from matbench_discovery.phonons import read_kappa_json
from scripts.evals import evaluate_models


def evaluate_one(model: Model) -> str | None:
    """Recompute one model's kappa metrics and update its YAML."""
    kappa_103_path = model.kappa_103_path
    if not kappa_103_path or not os.path.isfile(kappa_103_path):
        return "no kappa_103_path found"

    print(f"\nProcessing {model.label}...")
    df_ml = read_kappa_json(kappa_103_path)
    metrics = phonons.evaluate_kappa_predictions(df_ml)

    for key in ("srme", "sre", "srd"):
        print(f"\tκ_{key.upper()}={metrics[key]:.4f}")
    print(f"\tκ failure rate={metrics['failure_rate']:.2%}")
    print(f"\tImaginary-mode rate={metrics['imaginary_mode_rate']:.2%}")
    spectrum_w1 = metrics["spectrum_w1"]
    spectrum_w1_text = f"{spectrum_w1:.4f} THz" if spectrum_w1 is not None else "n/a"
    print(f"\tSpectrum W1={spectrum_w1_text}")

    phonons.write_metrics_to_yaml(model, metrics, kappa_103_path)
    print(f"\tUpdated {model.yaml_path}")
    return None


def main() -> int:
    """Evaluate kappa metrics, returning 0 if any model succeeds and 1 otherwise."""
    return evaluate_models("kappa", cli_args.models, evaluate_one)


if __name__ == "__main__":
    raise SystemExit(main())
