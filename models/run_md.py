"""Unified molecular dynamics benchmark runner for all MLIP models.

Every model registered in matbench_discovery.md_models runs through this one script.
Model dependency trees conflict, so each model resolves its own environment via
``uv run --with`` rather than sharing one. Typical cluster usage:

    # print the uv command that resolves <model>'s env and runs the full benchmark
    uv run models/run_md.py --print-cmd --model orb_v3

    # smoke-test a model end-to-end in seconds before committing to the long run
    uv run --with orb-models models/run_md.py --model orb_v3 --dry-run

    # full 20 ps NVT rollout + metrics over all CFPMD-26 systems (auto-downloaded)
    uv run --with orb-models models/run_md.py --model orb_v3 --write-yaml

    # one slurm array task per system: --systems "${SYSTEMS[$SLURM_ARRAY_TASK_ID]}"
"""

# /// script
# requires-python = ">=3.11"
# dependencies = ["ase>=3.27", "numpy", "pandas", "scipy", "tqdm", "matbench-discovery"]
#
# [tool.uv.sources]
# matbench-discovery = { path = "../", editable = true }
# ///

import argparse
import os
import shlex

from matbench_discovery import today
from matbench_discovery.md import run_md_benchmark
from matbench_discovery.md_models import MD_MODELS, load_calculator

module_dir = os.path.dirname(__file__)


def main() -> int:
    """Parse args and run (or dry-run) the MD benchmark for one model."""
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
        help="Smoke test: 1 system, a few MD steps, capped reference slice",
    )
    parser.add_argument("--out-dir", help="Defaults to models/<arch>/<today>-md-nvt")
    parser.add_argument("--ref-dir", help="Defaults to auto-downloaded CFPMD-26 set")
    parser.add_argument("--settings-csv")
    parser.add_argument("--ref-frame-interval-fs", type=float)
    parser.add_argument("--systems", nargs="*", help="Subset of system dir names")
    parser.add_argument("--n-steps", type=int, default=80_000)
    parser.add_argument("--time-step-fs", type=float, default=0.25)
    parser.add_argument("--record-interval", type=int, default=10)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--write-yaml", action="store_true", help="Write metrics to the model YAML"
    )
    args = parser.parse_args()

    if args.list_models:
        for key, spec in MD_MODELS.items():
            print(f"{key}: {', '.join(spec.deps) or '(core deps only)'}")
        return 0

    if not args.model:
        parser.error("--model is required (or pass --list-models)")
    if args.model not in MD_MODELS:
        parser.error(f"unknown --model {args.model!r}, see --list-models")

    if args.print_cmd:
        run_args = ["--model", args.model]
        if args.dry_run:
            run_args.append("--dry-run")
        print(
            shlex.join(MD_MODELS[args.model].uv_run_cmd("models/run_md.py", *run_args))
        )
        return 0

    from matbench_discovery.enums import Model

    # resolve the submission Model once (None for debug models like emt without an
    # enum entry); its YAML path gives the arch subdir (orb_v3 -> models/orb)
    try:
        model = Model.from_ref(args.model)
    except ValueError:
        model = None
    arch_dir = os.path.dirname(model.rel_path) if model else args.model
    out_dir = args.out_dir or f"{module_dir}/{arch_dir}/{today}-md-nvt"

    calculator = load_calculator(args.model)
    df_md = run_md_benchmark(
        calculator=calculator,
        model_key=args.model,
        out_dir=out_dir,
        ref_dir=args.ref_dir,
        settings_csv=args.settings_csv,
        ref_frame_interval_fs=args.ref_frame_interval_fs,
        systems=args.systems,
        n_steps=args.n_steps,
        time_step_fs=args.time_step_fs,
        record_interval=args.record_interval,
        seed=args.seed,
        dry_run=args.dry_run,
    )

    from matbench_discovery.metrics import md as md_metrics

    model_metrics = md_metrics.calc_md_metrics(df_md)
    print(f"\n{args.model} model-level MD metrics:")
    for key, val in model_metrics.items():
        print(f"  {key}: {val:.4f}" if isinstance(val, float) else f"  {key}: {val}")

    if args.write_yaml and not args.dry_run:
        if model is None:  # debug models like emt have no YAML to write to
            print(f"Skipping metrics.md write: {args.model!r} is not a Model")
            return 0
        # match the suffix run_md_benchmark adds for --systems subset runs so the
        # YAML pred_file points at the file actually written
        suffix = f"-{'-'.join(args.systems)}" if args.systems else ""
        csv_path = f"{out_dir}/{args.model}-md-metrics{suffix}.csv.gz"
        md_metrics.write_metrics_to_yaml(model, model_metrics, pred_file_path=csv_path)
        print(f"Wrote metrics.md to {model.yaml_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
