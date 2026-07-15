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
from functools import partial
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
from matbench_discovery.data import (
    artifact_date_from_prefix,
    artifact_filename,
    make_file_ref,
    update_yaml_file,
)
from matbench_discovery.discovery import (
    ARCHIVED_DISCOVERY_MODELS,
    COST_PROVENANCE_KEYS,
    DISCOVERY_PRED_COL,
    DiscoveryArtifacts,
    RelaxationSettings,
    dry_run_settings,
    load_wbm_atoms,
    merge_discovery_shards,
    run_discovery_shard,
    write_discovery_artifacts,
)
from matbench_discovery.enums import MbdKey, Model
from matbench_discovery.hpc import effective_shard_args as _effective_shard_args
from matbench_discovery.hpc import partition_material_ids
from matbench_discovery.hpc import slurm_shard_selection as _slurm_shard_selection
from matbench_discovery.runner_cli import dependency_run_args, resolve_sharded_prefix

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


def _resolve_output_paths(
    *,
    out_dir: str,
    settings: RelaxationSettings,
    dry_run: bool,
    shard_dir: str | None,
    pred_file: str | None,
    geo_opt_file: str | None,
) -> tuple[str, str, str]:
    """Reuse one prior shard directory and emit canonical final artifact names."""
    dry_suffix = "-dry-run" if dry_run else ""
    run_name = f"{today}-wbm-IS2RE-{settings.ase_optimizer}{dry_suffix}"
    default_prefix = f"{out_dir}/{run_name}"
    selected_shard_dir, artifact_prefix = resolve_sharded_prefix(
        default_prefix=default_prefix,
        prior_shard_pattern=(
            f"{out_dir}/*-wbm-IS2RE-{settings.ase_optimizer}{dry_suffix}-shards"
        ),
        task="discovery",
        shard_dir=shard_dir,
    )
    artifact_date = artifact_date_from_prefix(artifact_prefix, fallback=today)
    artifact_dir = f"{out_dir}/dry-run" if dry_run else out_dir
    discovery_file = artifact_filename(artifact_date, "discovery")
    geo_opt_file_name = artifact_filename(artifact_date, "geo_opt")
    return (
        selected_shard_dir,
        pred_file or f"{artifact_dir}/{discovery_file}",
        geo_opt_file or f"{artifact_dir}/{geo_opt_file_name}",
    )


def _print_cmd_args(args: argparse.Namespace, model_key: str) -> list[str]:
    """Forward meaningful CLI options into a dependency-isolated uv command."""
    return dependency_run_args(
        args,
        model_key,
        {
            "dtype": args.dtype,
            "device": args.device,
            "out-dir": args.out_dir,
            "shard-dir": args.shard_dir,
            "pred-file": args.pred_file,
            "geo-opt-file": args.geo_opt_file,
            "max-force": args.max_force,
            "max-steps": args.max_steps,
            "ase-optimizer": args.ase_optimizer,
            "cell-filter": args.cell_filter,
            "n-shards": args.n_shards,
            "shard-index": args.shard_index,
        },
        ("dry-run", "merge-shards", "write-yaml"),
    )


def _repo_relative_path(path: str) -> str:
    """Return a repository-relative path, preserving external absolute paths."""
    absolute_path = os.path.abspath(path)
    try:
        relative_path = os.path.relpath(absolute_path, ROOT).replace("\\", "/")
    except ValueError:  # Windows: path and repo root are on different drives
        return absolute_path
    return relative_path if not relative_path.startswith("../") else absolute_path


def _clear_legacy_and_update(
    section: dict[str, Any], *, updates: dict[str, Any]
) -> dict[str, Any]:
    """Drop stale flat artifact and cost keys, then apply regenerated metadata.

    Cost provenance is cleared before the merge so omitted fields (e.g. missing
    max_gpu_mem_gb) cannot linger from a prior YAML write.
    """
    legacy_keys = ("pred_file_url", "pred_file_artifact", "struct_col")
    for key in (*legacy_keys, *COST_PROVENANCE_KEYS):
        section.pop(key, None)
    section.update(updates)
    return section


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
    model_preds = artifacts.predictions[DISCOVERY_PRED_COL].copy()
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
    for task, path in (
        ("discovery", artifacts.pred_file_path),
        ("geo_opt", artifacts.geo_opt_file_path),
    ):
        artifact_data = {"pred_file": make_file_ref(_repo_relative_path(path))}
        if task == "discovery":
            artifact_data |= cost_data
        update_yaml_file(
            model.yaml_path,
            f"metrics.{task}",
            partial(_clear_legacy_and_update, updates=artifact_data),
        )
        print(f"Updated {task}.pred_file; re-upload if a prior Figshare URL existed")
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
        artifacts = write_discovery_artifacts(
            merged_run,
            pred_file_path=pred_file_path,
            geo_opt_file_path=geo_opt_file_path,
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
