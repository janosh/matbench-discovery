"""Unified diatomic pair-repulsion curve runner for all MLIP models.

Every model registered in matbench_discovery.calculators runs through this one script,
the same registry that powers ``models/run_md.py``. Model dependency trees conflict, so
each model resolves its own environment via ``uv run --with`` rather than sharing one.
Typical cluster usage:

    # print the uv command that resolves <model>'s env and runs the full sweep
    uv run models/run_diatomics.py --print-cmd --model orb_v3

    # smoke-test a model end-to-end in seconds before committing to the full sweep
    uv run --with orb-models models/run_diatomics.py --model orb_v3 --dry-run

    # full homonuclear sweep (H-U) + metrics written to the model YAML
    uv run --with orb-models models/run_diatomics.py --model orb_v3 --write-yaml
"""

# /// script
# requires-python = ">=3.11"
# dependencies = ["ase>=3.27", "numpy", "pandas", "scipy", "tqdm", "matbench-discovery"]
#
# [tool.uv.sources]
# matbench-discovery = { path = "../", editable = true }
# ///

import argparse
import gzip
import json
import os
import platform
import shlex
import time

import numpy as np

from matbench_discovery import today
from matbench_discovery.calculators import (
    CALCULATORS,
    load_calculator,
    resolve_calculator_key,
)
from matbench_discovery.diatomics import calc_diatomic_curve, homo_nuc

module_dir = os.path.dirname(__file__)
DIATOMIC_METRIC_EXCLUSIONS: dict[str, tuple[str, ...]] = {
    # These models deterministically blow up on He-He at short scored distances.
    "alphanet_v1_oam": ("He-He",),
    "sevennet_l3i5": ("He-He",),
}


def drop_metric_exclusions(
    model_key: str, metrics: dict[str, dict[str, float]]
) -> dict[str, dict[str, float]]:
    """Remove model-specific pathological curves before metric aggregation."""
    excluded = set(DIATOMIC_METRIC_EXCLUSIONS.get(model_key, ()))
    return {
        key: val
        for key, val in metrics.items()
        if key not in excluded and f"{key}-{key}" not in excluded
    }


def detect_hardware() -> str:
    """Human-readable name for the accelerator the run executed on.

    Tries torch, then JAX, then TensorFlow (MLIP backends expose the GPU differently),
    falling back to the CPU model when no GPU is visible.
    """
    try:
        import torch

        if torch.cuda.is_available():
            return torch.cuda.get_device_name(0)
    except ImportError:
        pass
    try:
        import jax

        if gpus := [dev for dev in jax.devices() if dev.platform == "gpu"]:
            return gpus[0].device_kind
    except ImportError:
        pass
    try:
        import tensorflow as tf

        if gpus := tf.config.list_physical_devices("GPU"):
            details = tf.config.experimental.get_device_details(gpus[0])
            return details.get("device_name", "GPU")
    except ImportError:
        pass
    return f"CPU ({platform.processor() or platform.machine()})"


def validate_distance_bounds(
    parser: argparse.ArgumentParser, min_dist: float, max_dist: float
) -> None:
    """Raise a clear CLI error for invalid geometric distance bounds."""
    if min_dist <= 0 or max_dist <= 0 or max_dist <= min_dist:
        parser.error("--min-dist and --max-dist must be positive with max > min")


def main() -> int:
    """Parse args and run (or dry-run) the diatomic curve sweep for one model."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", help="Model key, see --list-models")
    parser.add_argument(
        "--list-models", action="store_true", help="Print registered models and exit"
    )
    parser.add_argument(
        "--print-cmd",
        action="store_true",
        help="Print the 'uv run --with ...' command for --model and exit",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Smoke test: a few elements at a few distances",
    )
    parser.add_argument("--out-dir", help="Defaults to models/<arch>/<model>")
    parser.add_argument(
        "--max-z",
        type=int,
        default=92,
        help="Max atomic number for homonuclear pairs (default 92 = U)",
    )
    parser.add_argument(
        "--min-dist", type=float, default=0.1, help="Min separation in Å"
    )
    parser.add_argument(
        "--max-dist", type=float, default=6.0, help="Max separation in Å"
    )
    parser.add_argument(
        "--n-points", type=int, default=119, help="Number of log-spaced distances"
    )
    parser.add_argument(
        "--dtype",
        default="float64",
        choices=("float64", "float32"),
        help="Calculator float precision. default float64",
    )
    parser.add_argument(
        "--write-yaml",
        action="store_true",
        help="Compute + write metrics to model YAML",
    )
    args = parser.parse_args()

    if args.list_models:
        for key, spec in CALCULATORS.items():
            print(f"{key}: {', '.join(spec.deps) or '(core deps only)'}")
        return 0

    if not args.model:
        parser.error("--model is required (or pass --list-models)")
    try:
        args.model = resolve_calculator_key(args.model)
    except ValueError as exc:
        parser.error(f"{exc}, see --list-models")
    validate_distance_bounds(parser, args.min_dist, args.max_dist)
    if args.write_yaml and args.dry_run:
        parser.error("--write-yaml is incompatible with --dry-run")

    if args.print_cmd:
        run_args = [
            "--model",
            args.model,
            "--dtype",
            args.dtype,
            "--min-dist",
            str(args.min_dist),
            "--max-dist",
            str(args.max_dist),
            "--n-points",
            str(args.n_points),
            "--max-z",
            str(args.max_z),
        ]
        if args.out_dir:
            run_args += ["--out-dir", args.out_dir]
        if args.dry_run:
            run_args.append("--dry-run")
        if args.write_yaml:
            run_args.append("--write-yaml")
        cmd = CALCULATORS[args.model].uv_run_cmd("models/run_diatomics.py", *run_args)
        print(shlex.join(cmd))
        return 0

    from matbench_discovery.enums import Model

    # resolve the submission Model once (None for debug models like emt without an enum
    # entry); its YAML path stem gives the per-model dir (mace_mpa_0 -> mace/mace-mpa-0)
    try:
        model = Model.from_ref(args.model)
    except ValueError:
        model = None
    model_dir = os.path.splitext(model.rel_path)[0] if model else args.model
    out_dir = args.out_dir or f"{module_dir}/{model_dir}"

    max_z = 3 if args.dry_run else args.max_z
    n_points = 10 if args.dry_run else args.n_points
    # geometric spacing densifies the short-range repulsive wall, where unphysical
    # wiggles and discontinuities are most diagnostic (matches the MACE-MP convention)
    distances = np.geomspace(args.min_dist, args.max_dist, n_points)
    homo_pairs = [(z, z) for z in range(1, max_z + 1)]

    start_time = time.perf_counter()
    calculator = load_calculator(args.model, dtype=args.dtype)
    hardware = detect_hardware()  # after load: the backend's GPU is now initialized
    curves: dict[str, dict[str, list[float | list[list[float]]]]] = {}
    print(
        f"\nPredicting {len(homo_pairs)} diatomic curves for {args.model} on {hardware}"
    )
    calc_diatomic_curve(
        pairs=homo_pairs,
        calculator=calculator,
        model_name=args.model,
        distances=distances,
        results=curves,
    )
    run_time_sec = round(time.perf_counter() - start_time, 2)

    # drop curves with empty or non-finite energies/forces: NaN/Infinity are not valid
    # JSON (downstream JS parsers reject them), and these curves are skipped by the
    # metrics anyway, so the saved file matches what's actually scored
    def is_valid_curve(curve: dict[str, list[float | list[list[float]]]]) -> bool:
        """Return whether a curve has finite energy and force samples."""
        return bool(
            bool(curve.get("energies"))
            and bool(curve.get("forces"))
            and np.isfinite(curve["energies"]).all()
            and np.isfinite(curve["forces"]).all()
        )

    invalid_formulas = sorted(
        formula for formula, curve in curves.items() if not is_valid_curve(curve)
    )
    curves = {
        formula: curve for formula, curve in curves.items() if is_valid_curve(curve)
    }

    # flat on-disk schema read by DiatomicCurves.from_dict / the site's diatomics parser
    results = {homo_nuc: curves, "hetero-nuclear": {}, "distances": distances.tolist()}

    os.makedirs(out_dir, exist_ok=True)
    json_path = f"{out_dir}/{today}-diatomics.json.gz"
    with gzip.open(json_path, mode="wt") as file:
        # allow_nan=False guards against writing invalid JSON (NaN) if a non-finite
        # value slips through the filter above
        json.dump(results, file, allow_nan=False, default=lambda arr: arr.tolist())
    n_curves = len(curves)
    print(f"Wrote {n_curves}/{len(homo_pairs)} diatomic curves to {json_path}")
    print(f"Ran on {hardware} in {run_time_sec:.1f} s")

    if args.write_yaml:
        if model is None:  # debug models like emt have no YAML to write to
            print(f"Skipping metrics.diatomics write: {args.model!r} is not a Model")
            return 0
        from matbench_discovery.metrics import diatomics
        from matbench_discovery.metrics.diatomics import DiatomicCurves

        pred_curves = DiatomicCurves.from_dict(results)
        metrics = diatomics.calc_diatomic_metrics(
            ref_curves=diatomics.load_dft_reference_curves("PBE"),
            pred_curves=pred_curves,
            interpolate=200,
        )
        metrics = drop_metric_exclusions(args.model, metrics)
        excluded_formulas = sorted(
            {*DIATOMIC_METRIC_EXCLUSIONS.get(args.model, ()), *invalid_formulas}
        )
        run_metadata: dict[str, str | float | list[str]] = {
            "hardware": hardware,
            "run_time_sec": run_time_sec,
            "excluded_formulas": excluded_formulas,
        }
        mean_metrics = diatomics.write_metrics_to_yaml(
            model,
            metrics,
            pred_file_path=json_path,
            run_metadata=run_metadata,
        )
        print(f"\n{args.model} mean diatomic metrics:")
        for key, val in mean_metrics.items():
            print(f"  {key}: {val}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
