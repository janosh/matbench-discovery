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
import shlex
import time
from glob import glob

import numpy as np
from ase.data import chemical_symbols

from matbench_discovery import today
from matbench_discovery.calculators import (
    CALCULATORS,
    load_calculator,
    resolve_cli_calculator,
)
from matbench_discovery.data import artifact_filename
from matbench_discovery.diatomics import (
    CurveDict,
    DiatomicResults,
    calc_diatomic_curve,
    homo_nuc,
)
from matbench_discovery.hpc import detect_hardware, merge_run_metadata, peak_memory_gb
from matbench_discovery.metrics.diatomics import DIATOMIC_WALL_R_MIN_FACTOR, eval_window
from matbench_discovery.metrics.diatomics.exclusions import (
    drop_metric_exclusions,
    get_excluded_formula_reasons,
    is_non_mp_formula,
)
from matbench_discovery.runner_cli import add_common_runner_args, dependency_run_args

module_dir = os.path.dirname(__file__)


def trim_curve_to_finite(formula: str, curve: CurveDict) -> CurveDict | None:
    """Trim unscored non-finite points or reject a curve with scored ones.

    JSON cannot store NaN/Infinity. The scored range starts at 0.8x covalent radius;
    a non-finite value within it returns None so the curve is recorded as excluded.
    """
    energies = np.asarray(curve["energies"], dtype=float)
    forces = np.asarray(curve["forces"], dtype=float)
    if energies.size == 0 or forces.size == 0:
        return None
    distances = np.asarray(curve["distances"], dtype=float)
    finite = np.isfinite(energies) & np.isfinite(forces).all(axis=(1, 2))
    r_min, r_max = eval_window(
        formula,
        float(distances.max()),
        r_min_factor=DIATOMIC_WALL_R_MIN_FACTOR,
    )
    if not finite[(distances >= r_min - 1e-12) & (distances <= r_max)].all():
        return None
    if finite.all():
        return curve
    return {
        "distances": distances[finite].tolist(),
        "energies": energies[finite].tolist(),
        "forces": forces[finite].tolist(),
    }


def main() -> int:
    """Parse args and run (or dry-run) the diatomic curve sweep for one model."""
    parser = add_common_runner_args(
        argparse.ArgumentParser(description=__doc__),
        dry_run_help="Smoke test: a few elements at a few distances",
    )
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
        "--write-yaml",
        action="store_true",
        help="Compute + write metrics to model YAML",
    )
    parser.add_argument(
        "--merge-shards", action="store_true", help="Merge Slurm shards"
    )
    args = parser.parse_args()

    model_key = resolve_cli_calculator(parser, args.model, list_models=args.list_models)
    if model_key is None:
        return 0
    args.model = model_key
    if args.min_dist <= 0 or args.max_dist <= 0 or args.max_dist <= args.min_dist:
        parser.error("--min-dist and --max-dist must be positive with max > min")
    if args.max_z < 1:
        parser.error("--max-z must be a positive integer")
    if args.n_points < 2:
        parser.error("--n-points must be >= 2 (a curve needs at least 2 points)")
    if args.write_yaml and args.dry_run:
        parser.error("--write-yaml is incompatible with --dry-run")
    slurm_task_id = os.getenv("SLURM_ARRAY_TASK_ID")
    if args.write_yaml and slurm_task_id and not args.merge_shards:
        parser.error("--write-yaml is only supported by the --merge-shards task")

    if args.print_cmd:
        run_args = dependency_run_args(
            args,
            args.model,
            {
                "dtype": args.dtype,
                "min-dist": args.min_dist,
                "max-dist": args.max_dist,
                "n-points": args.n_points,
                "max-z": args.max_z,
                "out-dir": args.out_dir,
            },
            ("dry-run", "merge-shards", "write-yaml"),
        )
        print(
            shlex.join(
                CALCULATORS[args.model].uv_run_cmd("models/run_diatomics.py", *run_args)
            )
        )
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
    if slurm_task_id and not args.merge_shards:
        slurm_task_idx = int(slurm_task_id)
        if not 0 <= slurm_task_idx < max_z:
            parser.error(
                f"SLURM_ARRAY_TASK_ID must be in [0, {max_z - 1}] for {max_z=}"
            )
        z_values = [slurm_task_idx + 1]  # Slurm array indices are zero-based.
    else:
        z_values = range(1, max_z + 1)
    n_points = 10 if args.dry_run else args.n_points
    # geometric spacing densifies the short-range repulsive wall, where unphysical
    # wiggles and discontinuities are most diagnostic (matches the MACE-MP convention)
    distances = np.geomspace(args.min_dist, args.max_dist, n_points)
    shard_dir = f"{out_dir}/{today}-diatomics-shards"
    n_pairs = len(z_values)
    run_metadata: dict[str, str | float | dict[str, str]]
    curves: DiatomicResults = {}

    if args.merge_shards:
        shard_metadatas: list[dict[str, object]] = []
        # Normalize both sides because glob returns OS-native separators on Windows.
        expected_paths = [
            os.path.normpath(f"{shard_dir}/Z{z_value:03d}-diatomics.json.gz")
            for z_value in z_values
        ]
        shard_paths = sorted(
            map(os.path.normpath, glob(f"{shard_dir}/Z*-diatomics.json.gz"))
        )
        if shard_paths != expected_paths:
            parser.error(f"Expected shard files {expected_paths}, got {shard_paths}")
        for shard_path in shard_paths:
            with gzip.open(shard_path, mode="rt") as file:
                shard = json.load(file)
            shard_distances = shard["distances"]
            if shard_path == shard_paths[0]:
                distances = np.asarray(shard_distances)
            elif shard_distances != distances.tolist():
                parser.error(f"Inconsistent distances in {shard_path}")
            curves.update(shard.get(homo_nuc, {}))
            shard_metadatas.append(shard.get("run_metadata", {}))
        expected_formulas = {
            f"{chemical_symbols[z_value]}-{chemical_symbols[z_value]}"
            for z_value in z_values
        }
        missing_formulas = {
            formula
            for formula in expected_formulas - set(curves)
            # metrics skip non-MP elements, so models may omit their curves
            if not is_non_mp_formula(formula)
        }
        # deliberately no run-discovered invalid_formulas here (unlike the single-run
        # branch below): a shard's self-reported exclusions must not leak into the
        # merged output, so missing curves are only accepted once curated in the YAML
        exclusion_reasons = get_excluded_formula_reasons(args.model)
        if unexpected_missing := missing_formulas - set(exclusion_reasons):
            parser.error(f"Missing curves in shards: {sorted(unexpected_missing)}")
        run_metadata = {**merge_run_metadata(shard_metadatas)}
        json_path = f"{out_dir}/{artifact_filename(today, 'diatomics')}"
    else:
        start_time = time.perf_counter()
        calculator = load_calculator(args.model, dtype=args.dtype)
        hardware = detect_hardware()  # after load: the backend's GPU is now initialized
        print(f"\nPredicting {n_pairs} diatomic curves for {args.model} on {hardware}")
        calc_diatomic_curve(
            pairs=[(z_value, z_value) for z_value in z_values],
            calculator=calculator,
            model_name=args.model,
            distances=distances,
            results=curves,
        )
        run_time_sec = round(time.perf_counter() - start_time, 2)

        # trim out-of-window non-finite samples (unscored anyway, and NaN/Infinity are
        # not valid JSON); drop curves that are non-finite at scored separations
        invalid_formulas = []
        for formula in list(curves):
            trimmed = trim_curve_to_finite(formula, curves[formula])
            if trimmed is None:
                invalid_formulas.append(formula)
                del curves[formula]
            else:
                curves[formula] = trimmed

        exclusion_reasons = get_excluded_formula_reasons(args.model, invalid_formulas)
        # peak memory covers calculator load + full sweep (one process per run/shard)
        run_metadata = {
            "hardware": hardware,
            "run_time_sec": run_time_sec,
            **peak_memory_gb(),
        }
        if slurm_task_id:
            json_path = f"{shard_dir}/Z{z_values[0]:03d}-diatomics.json.gz"
        else:
            json_path = f"{out_dir}/{artifact_filename(today, 'diatomics')}"

    run_metadata["excluded_formula_reasons"] = exclusion_reasons
    # flat on-disk schema read by DiatomicCurves.from_dict / the site's diatomics parser
    results: dict[str, object] = {
        homo_nuc: curves,
        "hetero-nuclear": {},
        "distances": distances.tolist(),
    }
    # persist provenance (hardware, run time, exclusions) so shards can be aggregated
    # and merged/local files keep it
    results["run_metadata"] = run_metadata
    os.makedirs(os.path.dirname(json_path), exist_ok=True)
    with gzip.open(json_path, mode="wt") as file:
        # allow_nan=False guards against writing invalid JSON (NaN) if a non-finite
        # value slips through the filter above
        json.dump(results, file, allow_nan=False, default=lambda arr: arr.tolist())
    print(f"Wrote {len(curves)}/{n_pairs} diatomic curves to {json_path}")
    if not args.merge_shards:
        print(f"Ran on {hardware} in {run_time_sec:.1f} s")

    if args.write_yaml:
        if model is None:  # debug models like emt have no YAML to write to
            print(f"Skipping metrics.diatomics write: {args.model!r} is not a Model")
            return 0
        from matbench_discovery.metrics import diatomics

        metrics = diatomics.calc_diatomic_metrics(
            ref_curves=diatomics.load_dft_reference_curves("PBE"),
            pred_curves=diatomics.DiatomicCurves.from_dict(results),
            interpolate=200,
        )
        metrics = drop_metric_exclusions(args.model, metrics)
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
