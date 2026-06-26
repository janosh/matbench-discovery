"""Evaluate pre-computed MLIP MD trajectories against ab-initio references.

Computes per-system RDF, VDOS and (when stress is present) pressure metrics for
every model with trajectories in the MLIP directory and writes one per-system
metrics CSV per model. Energy/force RMSE requires running the model itself and is
handled by the unified runner models/run_md.py (which runs the model itself).

Expected data layout (CFPMD-26 convention):
    <ref-file>                               single reference HDF5 (group per system)
    <pred-dir>/<system>/nvt_<model>.extxyz   MLIP rollouts

Example:
    python scripts/evals/md_trajectories.py \
        --ref-file ~/data/cfpmd-26/cfpmd-26-aimd-reference.h5 \
        --pred-dir ~/data/cfpmd-26/mlip_trajectories
"""

import argparse
import os
import traceback
from glob import glob

import pandas as pd
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.md import (
    default_md_reference_path,
    read_reference_trajectory,
    read_trajectory,
)
from matbench_discovery.metrics import md as md_metrics
from matbench_discovery.trajectory import Trajectory

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--ref-file",
    help="Reference HDF5 (one group per system, with dt_fs/temperature_kelvin attrs). "
    "Defaults to the CFPMD-26 reference dataset, auto-downloaded from figshare.",
)
parser.add_argument("--pred-dir", required=True)
parser.add_argument(
    "--pred-frame-interval-fs",
    type=float,
    default=2.5,  # benchmark protocol: 0.25 fs timestep, frames every 10 steps
    help="Time between saved MLIP frames in fs",
)
parser.add_argument("--models", nargs="*", help="Subset of model keys to evaluate")
parser.add_argument("--out-dir", default=f"md-trajectory-metrics/{today}")
args = parser.parse_args()

ref_file = args.ref_file or default_md_reference_path()
os.makedirs(args.out_dir, exist_ok=True)

# discover (system, model) pairs from nvt_<model>.extxyz files in pred-dir
rows_by_model: dict[str, list[dict[str, float | str]]] = {}
pred_files = sorted(glob(f"{args.pred_dir}/*/nvt_*.extxyz*"))
# cache: each ref (trajectory, dt_fs) reused by ~15 models; keep at most one in memory
ref_cache: dict[str, tuple[Trajectory, float]] = {}

for pred_file in (pbar := tqdm(pred_files, desc="MD trajectory pairs")):
    system_name = os.path.basename(os.path.dirname(pred_file))
    model_key = os.path.basename(pred_file).removeprefix("nvt_").split(".extxyz")[0]
    if args.models and model_key not in args.models:
        continue
    pbar.set_postfix_str(f"{system_name} {model_key}")

    try:
        if system_name not in ref_cache:
            ref_cache.clear()  # keep at most one ref trajectory in memory
            ref_traj, ref_time_step_fs, _temperature = read_reference_trajectory(
                ref_file, system_name
            )
            ref_cache[system_name] = (ref_traj, ref_time_step_fs)
        ref_traj, ref_time_step_fs = ref_cache[system_name]
        system_metrics = md_metrics.evaluate_md_system(
            ref_traj,
            read_trajectory(pred_file),
            ref_time_step_fs=ref_time_step_fs,
            pred_time_step_fs=args.pred_frame_interval_fs,
        )
    except Exception as exc:  # noqa: BLE001 - record failure, keep evaluating
        print(f"\nFailed {system_name}/{model_key}: {exc!r}")
        traceback.print_exc()
        continue

    rows_by_model.setdefault(model_key, []).append(
        {"system": system_name, **system_metrics}
    )

if not rows_by_model:  # fail loudly instead of exiting 0 with no metrics written
    raise SystemExit(
        f"No trajectory pairs evaluated from {len(pred_files)} file(s) under "
        f"{args.pred_dir!r}: check the --pred-dir layout, the --models filter "
        f"({args.models}), and the per-pair errors printed above"
    )

for model_key, rows in rows_by_model.items():
    df_md = pd.DataFrame(rows).set_index("system")
    csv_path = f"{args.out_dir}/{today}-{model_key}-md-metrics.csv.gz"
    df_md.to_csv(csv_path)
    model_metrics = md_metrics.calc_md_metrics(df_md)
    metrics_str = ", ".join(f"{key}={val:.4g}" for key, val in model_metrics.items())
    print(f"{model_key}: {metrics_str}\n  saved to {csv_path}")
