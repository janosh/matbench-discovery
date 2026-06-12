"""Run NVT molecular dynamics with MACE and evaluate MD benchmark metrics.

Expects a reference directory with one subdirectory per system, each containing a
single ASE-readable trajectory of ab-initio MD frames. System directory names must
contain the simulation temperature as e.g. 'bulkAu_1500K'.

Example (auto-downloads the CFPMD-26 reference dataset on first use):
    python models/mace/test_mace_md.py
"""

import argparse
import json
import os
from glob import glob
from importlib.metadata import version

import ase.io
import pandas as pd
import torch
from mace.calculators import mace_mp
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.md import (
    default_md_reference_paths,
    extract_temperature,
    load_frame_intervals,
    read_trajectory,
    resolve_frame_interval,
    run_nvt_md,
)
from matbench_discovery.metrics import md as md_metrics

module_dir = os.path.dirname(__file__)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--ref-dir",
    help="Directory of per-system reference trajectories. Defaults to the CFPMD-26 "
    "reference dataset, auto-downloaded from figshare on first use.",
)
parser.add_argument("--model-name", default="mace-mp-0")
parser.add_argument(
    "--checkpoint", default="medium", help="MACE model name, path or URL"
)
interval_args = parser.add_mutually_exclusive_group()
interval_args.add_argument(
    "--ref-frame-interval-fs",
    type=float,
    help="Time between saved reference frames in fs (same for all systems)",
)
interval_args.add_argument(
    "--settings-csv",
    help="CSV with per-system reference frame dt. Defaults to the CFPMD-26 "
    "settings CSV when --ref-dir is also defaulted.",
)
parser.add_argument(
    "--systems",
    nargs="*",
    help="Subset of system dir names to evaluate (e.g. one per slurm array task)",
)
parser.add_argument("--n-steps", type=int, default=80_000)
parser.add_argument("--time-step-fs", type=float, default=0.25)
parser.add_argument("--record-interval", type=int, default=10)
parser.add_argument("--seed", type=int, default=0)
parser.add_argument(
    "--write-yaml", action="store_true", help="Write metrics to the model YAML file"
)
args = parser.parse_args()

if args.ref_dir is None:  # default to auto-downloaded CFPMD-26 reference dataset
    args.ref_dir, default_csv = default_md_reference_paths()
    args.settings_csv = args.settings_csv or default_csv
elif args.settings_csv is None and args.ref_frame_interval_fs is None:
    raise SystemExit("Pass --settings-csv or --ref-frame-interval-fs for custom data")

device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Loading MACE checkpoint {args.checkpoint!r} on {device=}")
calculator = mace_mp(
    model=args.checkpoint,
    device=device,
    default_dtype="float64",
    enable_cueq=device == "cuda",
)

out_dir = f"{module_dir}/{args.model_name}/{today}-md-nvt"
os.makedirs(out_dir, exist_ok=True)

run_params = dict(
    **vars(args),
    device=device,
    versions={dep: version(dep) for dep in ("ase", "mace-torch", "torch")},
)
with open(f"{out_dir}/run_params.json", mode="w") as file:
    json.dump(run_params, file, indent=4)

frame_intervals = load_frame_intervals(args.settings_csv) if args.settings_csv else {}
system_dirs = sorted(
    entry.path
    for entry in os.scandir(args.ref_dir)
    if entry.is_dir() and (not args.systems or entry.name in args.systems)
)
per_system_rows: list[dict[str, float | str]] = []

for system_dir in (pbar := tqdm(system_dirs, desc="MD systems")):
    system_name = os.path.basename(system_dir)
    pbar.set_postfix_str(system_name)

    # match traj.extxyz as well as ASE-transparent .gz/.xz/.bz2 compressed variants
    traj_files = sorted(glob(f"{system_dir}/traj.*xyz*"))
    if len(traj_files) != 1:
        raise ValueError(
            f"Expected 1 trajectory file in {system_dir}, got {traj_files}"
        )
    temperature_kelvin = extract_temperature(system_name)
    if temperature_kelvin is None:
        raise ValueError(f"Could not parse temperature from {system_name!r}")
    ref_frame_interval_fs = args.ref_frame_interval_fs or resolve_frame_interval(
        system_name, frame_intervals
    )
    if ref_frame_interval_fs is None:
        raise ValueError(f"No frame interval for {system_name!r} in settings CSV")

    ref_trajectory = read_trajectory(traj_files[0])
    pred_traj_path = f"{out_dir}/{system_name}-nvt-{args.model_name}.extxyz"
    if os.path.isfile(pred_traj_path):
        pred_trajectory = read_trajectory(pred_traj_path)
    else:
        pred_trajectory = run_nvt_md(
            ref_trajectory[0],
            calculator,
            temperature_kelvin=temperature_kelvin,
            n_steps=args.n_steps,
            time_step_fs=args.time_step_fs,
            record_interval=args.record_interval,
            seed=args.seed,
        )
        ase.io.write(pred_traj_path, pred_trajectory)

    system_metrics = md_metrics.evaluate_md_system(
        ref_trajectory,
        pred_trajectory,
        ref_time_step_fs=ref_frame_interval_fs,
        pred_time_step_fs=args.time_step_fs * args.record_interval,
        calculator=calculator,
    )
    per_system_rows.append(
        {"system": system_name, "temperature_kelvin": temperature_kelvin}
        | system_metrics
    )

df_md = pd.DataFrame(per_system_rows).set_index("system")
# suffix avoids collisions when parallel single-system jobs share out_dir
suffix = f"-{'-'.join(args.systems)}" if args.systems else ""
csv_path = f"{out_dir}/{today}-{args.model_name}-md-metrics{suffix}.csv.gz"
df_md.to_csv(csv_path)
print(f"Per-system metrics saved to {csv_path!r}")

model_metrics = md_metrics.calc_md_metrics(df_md)
print(json.dumps(model_metrics, indent=2))

if args.write_yaml:
    from matbench_discovery.enums import Model

    model = Model.from_ref(args.model_name)
    md_metrics.write_metrics_to_yaml(model, model_metrics, pred_file_path=csv_path)
    print(f"Updated {model.yaml_path}")
