"""Unified dependency-isolated WBM discovery runner for ASE calculators.

Examples:
    uv run models/run_discovery.py --list-models
    uv run models/run_discovery.py --print-cmd --model mace_mp_0 --dry-run
    uv run --with mace-torch models/run_discovery.py --model mace_mp_0
    uv run models/run_discovery.py --model mace_mp_0 --merge-shards --write-yaml
"""

# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "ase>=3.28",
#   "numpy",
#   "pandas",
#   "pymatgen>=2026.5.4",
#   "tqdm",
#   "matbench-discovery",
# ]
#
# [tool.uv.sources]
# matbench-discovery = { path = "../", editable = true }
# ///

from __future__ import annotations

import argparse
import glob
import os
import shlex
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

import pandas as pd

from matbench_discovery import ROOT, today
from matbench_discovery.calculators import (
    CALCULATORS,
    load_calculator,
    resolve_cli_calculator,
)
from matbench_discovery.data import update_yaml_file
from matbench_discovery.discovery import (
    ARCHIVED_DISCOVERY_MODELS,
    COST_PROVENANCE_KEYS,
    DiscoveryArtifacts,
    RelaxationSettings,
    dry_run_settings,
    load_wbm_atoms,
    merge_discovery_shards,
    partition_material_ids,
    run_discovery_shard,
    write_discovery_artifacts,
)
from matbench_discovery.enums import MbdKey, Model

module_dir = os.path.dirname(__file__)


def build_parser() -> argparse.ArgumentParser:
    """Build the discovery runner command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", help="Calculator key, model key, or model label")
    parser.add_argument(
        "--list-models", action="store_true", help="Print registered calculators"
    )
    parser.add_argument(
        "--print-cmd",
        action="store_true",
        help="Print the dependency-isolated uv command and exit",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Relax four WBM structures for at most two optimizer steps",
    )
    parser.add_argument(
        "--merge-shards",
        action="store_true",
        help="Strictly merge complete shards and write final artifacts",
    )
    parser.add_argument(
        "--write-yaml",
        action="store_true",
        help="On a complete merge, update artifact paths and discovery metrics",
    )
    parser.add_argument("--out-dir", help="Defaults to models/<arch>/<model>")
    parser.add_argument("--shard-dir", help="Override the resumable shard directory")
    parser.add_argument("--pred-file", help="Override the final prediction CSV path")
    parser.add_argument(
        "--geo-opt-file", help="Override the final relaxed-structure JSONL path"
    )
    parser.add_argument("--dtype", choices=("float64", "float32"), default="float64")
    parser.add_argument("--device", choices=("cpu", "cuda"))
    parser.add_argument("--max-force", type=float)
    parser.add_argument("--max-steps", type=int)
    parser.add_argument("--ase-optimizer")
    parser.add_argument(
        "--cell-filter",
        help="ASE filter class name, or 'none' for fixed-cell relaxation",
    )
    parser.add_argument(
        "--n-shards",
        type=int,
        help="Number of atom-balanced shards; defaults to Slurm task count or 1. "
        "When rerunning a subset of shards, pass the original value so task IDs "
        "keep their original shard indices",
    )
    parser.add_argument(
        "--shard-index",
        type=int,
        help="Zero-based shard index; defaults to normalized Slurm array task ID",
    )
    return parser


def _artifact_columns(model: Model | None) -> tuple[str | None, str | None]:
    """Preserve established prediction and structure columns from model YAML."""
    if model is None:
        return None, None

    def _metric_column(task: str, key: str) -> str | None:
        """Return one configured artifact column, if present."""
        metrics = model.metrics.get(task)
        value = metrics.get(key) if isinstance(metrics, dict) else None
        return str(value) if value else None

    return _metric_column("discovery", "pred_col"), _metric_column(
        "geo_opt", "struct_col"
    )


def _slurm_shard_selection(
    n_shards_arg: int | None, shard_index_arg: int | None
) -> tuple[int, int]:
    """Resolve zero-based shard selection from CLI flags or Slurm environment."""
    n_shards = (
        n_shards_arg
        if n_shards_arg is not None
        else int(os.getenv("SLURM_ARRAY_TASK_COUNT") or 1)
    )
    if n_shards < 1:
        raise ValueError(f"n_shards must be positive, got {n_shards}")

    if shard_index_arg is not None:
        shard_index = shard_index_arg
    elif slurm_task_id := os.getenv("SLURM_ARRAY_TASK_ID"):
        if n_shards_arg is not None:
            # Explicit shard counts make partial reruns such as --array=3,7 map back
            # to their original zero-based shard indices.
            slurm_task_max = os.getenv("SLURM_ARRAY_TASK_MAX")
            if slurm_task_max is not None and int(slurm_task_max) >= n_shards:
                # fail every task of a 1-based --array=1-N immediately instead of
                # running N-1 shards and silently never producing shard 0
                raise ValueError(
                    f"SLURM_ARRAY_TASK_MAX={slurm_task_max} exceeds the last shard "
                    f"index {n_shards - 1}; with explicit --n-shards, array task IDs "
                    "are treated as zero-based original shard indices"
                )
            shard_index = int(slurm_task_id)
        else:
            slurm_task_min = int(os.getenv("SLURM_ARRAY_TASK_MIN", "0"))
            slurm_task_max_id = int(
                os.getenv("SLURM_ARRAY_TASK_MAX", str(slurm_task_min + n_shards - 1))
            )
            if slurm_task_max_id - slurm_task_min + 1 != n_shards:
                raise ValueError(
                    "Partial Slurm arrays require explicit --n-shards so task IDs "
                    "retain their original shard indices (any rerun of a shard "
                    "subset, contiguous or not, must pass the original --n-shards)"
                )
            shard_index = int(slurm_task_id) - slurm_task_min
    elif n_shards == 1:
        shard_index = 0
    else:
        raise ValueError("--shard-index is required outside a Slurm array")

    if not 0 <= shard_index < n_shards:
        raise ValueError(
            f"shard_index must be in [0, {n_shards - 1}], got {shard_index}"
        )
    return n_shards, shard_index


def _effective_shard_args(
    n_shards: int | None,
    shard_index: int | None,
    *,
    dry_run: bool,
) -> tuple[int | None, int | None, bool]:
    """Run one dry-run shard on only the first task of an implicit Slurm array."""
    if (
        not dry_run
        or n_shards is not None
        or shard_index is not None
        or not os.getenv("SLURM_ARRAY_TASK_COUNT")
    ):
        return n_shards, shard_index, False
    slurm_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "0"))
    slurm_task_min = int(os.getenv("SLURM_ARRAY_TASK_MIN", "0"))
    return 1, 0, slurm_task_id != slurm_task_min


def _resolve_output_paths(
    *,
    out_dir: str,
    settings: RelaxationSettings,
    dry_run: bool,
    shard_dir: str | None,
    pred_file: str | None,
    geo_opt_file: str | None,
) -> tuple[str, str, str]:
    """Reuse one prior shard directory and keep final artifact dates aligned."""
    dry_suffix = "-dry-run" if dry_run else ""
    run_name = f"{today}-wbm-IS2RE-{settings.ase_optimizer}{dry_suffix}"
    default_prefix = f"{out_dir}/{run_name}"
    selected_shard_dir = shard_dir or f"{default_prefix}-shards"
    if shard_dir is None and not os.path.isdir(selected_shard_dir):
        pattern = f"{out_dir}/*-wbm-IS2RE-{settings.ase_optimizer}{dry_suffix}-shards"
        candidates = sorted(path for path in glob.glob(pattern) if os.path.isdir(path))
        if len(candidates) > 1:
            raise ValueError(
                f"Multiple discovery shard directories found: {candidates}. "
                "Select one with --shard-dir."
            )
        if candidates:
            selected_shard_dir = candidates[0]

    artifact_prefix = (
        selected_shard_dir.removesuffix("-shards")
        if selected_shard_dir.endswith("-shards")
        else default_prefix
    )
    return (
        selected_shard_dir,
        pred_file or f"{artifact_prefix}.csv.gz",
        geo_opt_file or f"{artifact_prefix}.jsonl.gz",
    )


def _print_cmd_args(args: argparse.Namespace, model_key: str) -> list[str]:
    """Forward meaningful CLI options into a dependency-isolated uv command."""
    run_args = ["--model", model_key, "--dtype", args.dtype]
    value_options = (
        ("device", args.device),
        ("out-dir", args.out_dir),
        ("shard-dir", args.shard_dir),
        ("pred-file", args.pred_file),
        ("geo-opt-file", args.geo_opt_file),
        ("max-force", args.max_force),
        ("max-steps", args.max_steps),
        ("ase-optimizer", args.ase_optimizer),
        ("cell-filter", args.cell_filter),
        ("n-shards", args.n_shards),
        ("shard-index", args.shard_index),
    )
    for option, value in value_options:
        if value is not None:
            run_args.extend([f"--{option}", str(value)])
    run_args.extend(
        f"--{flag}"
        for flag in ("dry-run", "merge-shards", "write-yaml")
        if getattr(args, flag.replace("-", "_"))
    )
    return run_args


def _repo_relative_path(path: str) -> str:
    """Return a repository-relative path, preserving external absolute paths."""
    absolute_path = os.path.abspath(path)
    try:
        relative_path = os.path.relpath(absolute_path, ROOT).replace("\\", "/")
    except ValueError:  # Windows: path and repo root are on different drives
        return absolute_path
    return relative_path if not relative_path.startswith("../") else absolute_path


def _artifact_yaml_data(
    model: Model,
    task: str,
    path: str,
    column_key: str,
    column: str,
) -> dict[str, str | float | None]:
    """Build artifact metadata, invalidating any previously uploaded URL.

    The local artifact was just regenerated, so an existing URL points at stale
    content even when the file path is unchanged (e.g. a redone shard re-merged
    on the same day).
    """
    data: dict[str, str | float | None] = {
        "pred_file": _repo_relative_path(path),
        column_key: column,
    }
    existing = model.metrics.get(task)
    if isinstance(existing, dict) and existing.get("pred_file_url"):
        data["pred_file_url"] = None
    return data


def _write_yaml_results(
    model: Model,
    artifacts: DiscoveryArtifacts,
    run_metadata: Mapping[str, Any] | None = None,
) -> None:
    """Update artifact paths, cost provenance, and all three discovery subsets."""
    from matbench_discovery.data import MAX_E_FORM_ERROR_THRESHOLD, df_wbm
    from matbench_discovery.metrics import discovery as discovery_metrics

    # mirror load_df_wbm_with_preds's outlier masking and .round(3) convention so
    # metrics written here match a later scripts/evals/discovery.py recompute
    model_preds = artifacts.predictions[artifacts.pred_col].copy()
    bad_mask = abs(model_preds - df_wbm[MbdKey.e_form_dft]) > MAX_E_FORM_ERROR_THRESHOLD
    model_preds.loc[bad_mask] = pd.NA
    metric_reference = df_wbm.round(3)
    model_preds = pd.to_numeric(
        model_preds.round(3).reindex(metric_reference.index), errors="coerce"
    )
    subset_indices = discovery_metrics.discovery_subset_indices(
        metric_reference, model_preds
    )
    metrics_by_subset = discovery_metrics.calc_discovery_metrics(
        metric_reference,
        model_preds,
        subset_indices=subset_indices,
        uniq_proto_prevalence=discovery_metrics.wbm_uniq_proto_prevalence(),
    )

    # rollout cost provenance, mirroring the MD and diatomics runners (fields are
    # dropped upstream unless every shard reported them on uniform hardware)
    cost_data = {
        key: value
        for key in COST_PROVENANCE_KEYS
        if (value := (run_metadata or {}).get(key)) is not None
    }
    for task, path, column_key, column in (
        ("discovery", artifacts.pred_file_path, "pred_col", artifacts.pred_col),
        ("geo_opt", artifacts.geo_opt_file_path, "struct_col", artifacts.struct_col),
    ):
        artifact_data = _artifact_yaml_data(model, task, path, column_key, column)
        if task == "discovery":
            artifact_data |= cost_data
        update_yaml_file(model.yaml_path, f"metrics.{task}", artifact_data)
        if "pred_file_url" in artifact_data:
            print(f"Cleared stale {task}.pred_file_url; upload the new artifact")
    discovery_metrics.write_all_metrics_to_yaml(
        model,
        metrics_by_subset,
        metric_reference,
        model_preds,
        subset_indices=subset_indices,
    )


def main(raw_args: Sequence[str] | None = None) -> int:
    """Run, resume, or merge the unified WBM discovery benchmark."""
    parser = build_parser()
    args = parser.parse_args(raw_args)

    model_key = resolve_cli_calculator(
        parser,
        args.model,
        list_models=args.list_models,
        archived_reasons=ARCHIVED_DISCOVERY_MODELS,
        task="discovery",
    )
    if model_key is None:
        return 0
    if args.write_yaml and not args.merge_shards:
        parser.error("--write-yaml is only supported with --merge-shards")
    if args.write_yaml and args.dry_run:
        parser.error("--write-yaml is incompatible with --dry-run")
    if args.merge_shards and args.shard_index is not None:
        parser.error("--shard-index is incompatible with --merge-shards")

    if args.print_cmd:
        command = CALCULATORS[model_key].uv_run_cmd(
            "models/run_discovery.py", *_print_cmd_args(args, model_key)
        )
        print(shlex.join(command))
        return 0

    try:
        settings = RelaxationSettings.from_model(
            model_key,
            max_force=args.max_force,
            max_steps=args.max_steps,
            ase_optimizer=args.ase_optimizer,
            cell_filter=args.cell_filter,
        )
        if args.dry_run:
            settings = dry_run_settings(settings)
    except (TypeError, ValueError) as exc:
        parser.error(str(exc))

    try:
        model = Model.from_ref(model_key)
    except ValueError:
        model = None
    model_dir = os.path.splitext(model.rel_path)[0] if model else model_key
    out_dir = args.out_dir or f"{module_dir}/{model_dir}"
    try:
        shard_dir, pred_file_path, geo_opt_file_path = _resolve_output_paths(
            out_dir=out_dir,
            settings=settings,
            dry_run=args.dry_run,
            shard_dir=args.shard_dir,
            pred_file=args.pred_file,
            geo_opt_file=args.geo_opt_file,
        )
    except ValueError as exc:
        parser.error(str(exc))

    if args.merge_shards:
        shard_suffix = f"{args.n_shards:03d}" if args.n_shards is not None else "*"
        shard_paths = sorted(glob.glob(f"{shard_dir}/shard-*-of-{shard_suffix}.jsonl"))
        if args.dry_run:
            expected_material_ids = list(load_wbm_atoms(model_key, dry_run=True))
        else:
            from matbench_discovery.data import df_wbm

            expected_material_ids = list(map(str, df_wbm.index))
        try:
            merged_run = merge_discovery_shards(
                shard_paths,
                model_key=model_key,
                expected_material_ids=expected_material_ids,
            )
        except (OSError, TypeError, ValueError) as exc:
            parser.error(str(exc))
        if args.n_shards is not None and merged_run.n_shards != args.n_shards:
            parser.error(
                f"Merged {merged_run.n_shards} shards, expected {args.n_shards}"
            )
        pred_col, struct_col = _artifact_columns(model)
        artifacts = write_discovery_artifacts(
            merged_run,
            pred_file_path=pred_file_path,
            geo_opt_file_path=geo_opt_file_path,
            pred_col=pred_col,
            struct_col=struct_col,
        )
        print(
            f"Wrote {artifacts.n_success:,} predictions "
            f"({artifacts.n_failed:,} failures) to {artifacts.pred_file_path}"
        )
        print(f"Wrote relaxed structures to {artifacts.geo_opt_file_path}")
        if args.write_yaml:
            if model is None:
                print(f"Skipping YAML write: {model_key!r} is not a Model")
            else:
                _write_yaml_results(model, artifacts, merged_run.run_metadata)
                print(f"Updated discovery and geo-opt metadata in {model.yaml_path}")
        return 0

    try:
        n_shards_arg, shard_index_arg, skip_task = _effective_shard_args(
            args.n_shards, args.shard_index, dry_run=args.dry_run
        )
        if skip_task:
            print("Skipping redundant dry run outside the first Slurm array task")
            return 0
        n_shards, shard_index = _slurm_shard_selection(n_shards_arg, shard_index_arg)
        atoms_by_id = load_wbm_atoms(model_key, dry_run=args.dry_run)
        shard_assignments = partition_material_ids(atoms_by_id, n_shards)
    except ValueError as exc:
        parser.error(str(exc))
    dataset_material_ids = list(atoms_by_id)
    assigned_material_ids = shard_assignments[shard_index]
    atoms_by_id = {
        material_id: atoms_by_id[material_id] for material_id in assigned_material_ids
    }
    shard_path = f"{shard_dir}/shard-{shard_index:03d}-of-{n_shards:03d}.jsonl"
    calculator = load_calculator(model_key, device=args.device, dtype=args.dtype)
    shard = run_discovery_shard(
        calculator=calculator,
        model_key=model_key,
        atoms_by_id=atoms_by_id,
        assigned_material_ids=assigned_material_ids,
        shard_path=shard_path,
        shard_index=shard_index,
        n_shards=n_shards,
        settings=settings,
        dataset_material_ids=dataset_material_ids,
    )
    n_failures = sum(record.error is not None for record in shard.records)
    print(
        f"Wrote shard {shard_index + 1}/{n_shards} with "
        f"{len(shard.records):,} records ({n_failures:,} failures) to {shard_path}"
    )
    return int(n_failures > 0)


if __name__ == "__main__":
    raise SystemExit(main())
