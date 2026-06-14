"""Evaluate pre-computed MLIP MD trajectories against ab-initio references.

Computes per-system RDF, VDOS and (when stress is present) pressure metrics for
every model with trajectories in the MLIP directory and writes one per-system
metrics CSV per model. Energy/force RMSE requires running the model itself and is
handled by the unified runner models/run_md.py (which runs the model itself).

Expected data layout (CFPMD-26 convention):
    <ref-dir>/<system>/traj.extxyz      reference AIMD trajectories
    <pred-dir>/<system>/nvt_<model>.extxyz   MLIP rollouts

Example:
    python scripts/evals/md_trajectories.py \
        --ref-dir ~/data/cfpmd-26/reference_AIMD_trajectories \
        --pred-dir ~/data/cfpmd-26/mlip_trajectories \
        --settings-csv ~/data/cfpmd-26/reference_AIMD_timestep_and_stride.csv
"""

import argparse
import os
import traceback
from glob import glob

import pandas as pd
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.md import (
    default_md_reference_paths,
    load_frame_intervals,
    read_trajectory,
    resolve_frame_interval,
)
from matbench_discovery.metrics import md as md_metrics

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--ref-dir",
    help="Directory of per-system reference trajectories. Defaults to the CFPMD-26 "
    "reference dataset, auto-downloaded from figshare on first use.",
)
parser.add_argument("--pred-dir", required=True)
parser.add_argument("--settings-csv", help="CSV with per-system reference frame dt")
parser.add_argument(
    "--pred-frame-interval-fs",
    type=float,
    default=2.5,  # benchmark protocol: 0.25 fs timestep, frames every 10 steps
    help="Time between saved MLIP frames in fs",
)
parser.add_argument("--models", nargs="*", help="Subset of model keys to evaluate")
parser.add_argument("--out-dir", default=f"md-trajectory-metrics/{today}")
args = parser.parse_args()

if args.ref_dir is None:  # default to auto-downloaded CFPMD-26 reference dataset
    args.ref_dir, default_csv = default_md_reference_paths()
    args.settings_csv = args.settings_csv or default_csv
elif args.settings_csv is None:
    raise SystemExit("Pass --settings-csv when using a custom --ref-dir")

frame_intervals = load_frame_intervals(args.settings_csv)
os.makedirs(args.out_dir, exist_ok=True)

# discover (system, model) pairs from nvt_<model>.extxyz files in pred-dir
rows_by_model: dict[str, list[dict[str, float | str]]] = {}
pred_files = sorted(glob(f"{args.pred_dir}/*/nvt_*.extxyz*"))
ref_trajectories: dict[str, list] = {}  # cache: each ref is reused by ~15 models

for pred_file in (pbar := tqdm(pred_files, desc="MD trajectory pairs")):
    system_name = os.path.basename(os.path.dirname(pred_file))
    model_key = os.path.basename(pred_file).removeprefix("nvt_").split(".extxyz")[0]
    if args.models and model_key not in args.models:
        continue
    pbar.set_postfix_str(f"{system_name} {model_key}")

    ref_dt_fs = resolve_frame_interval(system_name, frame_intervals)
    if ref_dt_fs is None:
        raise ValueError(f"No frame interval for {system_name!r} in settings CSV")

    try:
        if system_name not in ref_trajectories:
            # match traj.extxyz and ASE-transparent .gz/.xz/.bz2 compressed variants
            ref_files = glob(f"{args.ref_dir}/{system_name}/traj.*xyz*")
            if len(ref_files) != 1:
                raise ValueError(
                    f"Expected 1 reference trajectory for {system_name}, "
                    f"got {ref_files}"
                )
            ref_trajectories.clear()  # keep at most one ref trajectory in memory
            ref_trajectories[system_name] = read_trajectory(ref_files[0])
        system_metrics = md_metrics.evaluate_md_system(
            ref_trajectories[system_name],
            read_trajectory(pred_file),
            ref_time_step_fs=ref_dt_fs,
            pred_time_step_fs=args.pred_frame_interval_fs,
        )
    except Exception as exc:  # noqa: BLE001 - record failure, keep evaluating
        print(f"\nFailed {system_name}/{model_key}: {exc!r}")
        traceback.print_exc()
        continue

    rows_by_model.setdefault(model_key, []).append(
        {"system": system_name, **system_metrics}
    )

for model_key, rows in rows_by_model.items():
    df_md = pd.DataFrame(rows).set_index("system")
    csv_path = f"{args.out_dir}/{today}-{model_key}-md-metrics.csv.gz"
    df_md.to_csv(csv_path)
    model_metrics = md_metrics.calc_md_metrics(df_md)
    metrics_str = ", ".join(f"{key}={val:.4g}" for key, val in model_metrics.items())
    print(f"{model_key}: {metrics_str}\n  saved to {csv_path}")
