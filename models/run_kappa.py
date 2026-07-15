"""Unified dependency-isolated PhononDB kappa runner for registered MLIPs.

Examples:
    uv run models/run_kappa.py --list-models
    uv run models/run_kappa.py --print-cmd --model mace_mp_0 --dry-run
    uv run models/run_kappa.py --model mace_mp_0 --n-shards 8 --shard-index 0
    uv run models/run_kappa.py --model mace_mp_0 --merge-shards --write-yaml
"""

# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "matbench-discovery[phonons,symmetry]",
# ]
#
# [tool.uv.sources]
# matbench-discovery = { path = "../", editable = true }
# ///

from __future__ import annotations

import argparse
import os
import shlex
from typing import TYPE_CHECKING

from matbench_discovery import today
from matbench_discovery.calculators import (
    CALCULATORS,
    load_calculator,
    resolve_calculator_key,
    resolve_checkpoint,
    resolve_device,
)
from matbench_discovery.data import artifact_date_from_prefix, artifact_filename
from matbench_discovery.enums import DataFiles, Model
from matbench_discovery.hpc import effective_shard_args, slurm_shard_selection
from matbench_discovery.phonons.pipeline import (
    KappaArtifacts,
    KappaSettings,
    MergedKappaRun,
    checkpoint_digest,
    file_sha256,
    load_phonondb_atoms,
    merge_kappa_shards,
    run_kappa_shard,
    write_kappa_artifacts,
)
from matbench_discovery.runner_cli import dependency_run_args, resolve_sharded_prefix

if TYPE_CHECKING:
    from collections.abc import Sequence


def build_parser() -> argparse.ArgumentParser:
    """Build the unified kappa command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", help="Calculator key, model key, or model label")
    parser.add_argument(
        "--list-models",
        action="store_true",
        help="Print calculators with verified hyperparams.kappa settings",
    )
    parser.add_argument(
        "--print-cmd",
        action="store_true",
        help="Print the dependency-isolated uv command and exit",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Run a bounded one-structure end-to-end smoke test",
    )
    parser.add_argument(
        "--dry-run-size",
        type=int,
        default=1,
        help="Number of structures retained by --dry-run",
    )
    parser.add_argument(
        "--merge-shards",
        action="store_true",
        help="Strictly merge complete records and write final artifacts",
    )
    parser.add_argument(
        "--write-yaml",
        action="store_true",
        help="On a complete 103-structure merge, update metrics and provenance",
    )
    parser.add_argument(
        "--retry-failures",
        action="store_true",
        help="Recompute existing records that contain a persisted error",
    )
    parser.add_argument("--out-dir", help="Defaults to models/<arch>/<model>")
    parser.add_argument("--shard-dir", help="Override the resumable run directory")
    parser.add_argument(
        "--dataset",
        help="Override the canonical PhononDB 103-structure extxyz path",
    )
    parser.add_argument(
        "--checkpoint",
        help="Override the model checkpoint/artifact for supported calculators",
    )
    parser.add_argument("--dtype", choices=("float64", "float32"), default="float64")
    parser.add_argument("--device", choices=("cpu", "cuda"))
    parser.add_argument(
        "--n-shards",
        type=int,
        help="Number of atom-balanced shards; defaults to Slurm task count or 1",
    )
    parser.add_argument(
        "--shard-index",
        type=int,
        help="Zero-based shard index; defaults to normalized Slurm array task ID",
    )
    return parser


def kappa_model_keys() -> tuple[str, ...]:
    """Return registered calculators with versioned YAML-backed kappa settings."""
    model_keys = []
    for model_key in CALCULATORS:
        try:
            KappaSettings.from_model(model_key)
        except (TypeError, ValueError):
            continue
        model_keys.append(model_key)
    return tuple(model_keys)


def resolve_output_paths(
    *,
    out_dir: str,
    dry_run: bool,
    shard_dir: str | None,
) -> tuple[str, str]:
    """Reuse one prior shard directory and emit a canonical prediction artifact."""
    dry_suffix = "-dry-run" if dry_run else ""
    default_stem = f"{out_dir}/{today}-phonondb-kappa-103{dry_suffix}"
    selected_shard_dir, artifact_stem = resolve_sharded_prefix(
        default_prefix=default_stem,
        prior_shard_pattern=f"{out_dir}/*-phonondb-kappa-103{dry_suffix}-shards",
        task="kappa",
        shard_dir=shard_dir,
    )
    artifact_date = artifact_date_from_prefix(artifact_stem, fallback=today)
    artifact_dir = f"{out_dir}/dry-run" if dry_run else out_dir
    artifact_name = artifact_filename(artifact_date, "phonons_kappa_103")
    return selected_shard_dir, f"{artifact_dir}/{artifact_name}"


def resolved_artifact_identifier(model: Model, checkpoint: str | None) -> str | None:
    """Return an actual checkpoint path when managed, else a stable model ID."""
    checkpoint_path = resolve_checkpoint(model.name, checkpoint)
    if checkpoint_path:
        return checkpoint_path
    return model.metadata["checkpoint_url"]


def print_cmd_args(args: argparse.Namespace, model_key: str) -> list[str]:
    """Forward meaningful options into a dependency-isolated uv command."""
    return dependency_run_args(
        args,
        model_key,
        {
            "dtype": args.dtype,
            "device": args.device,
            "out-dir": args.out_dir,
            "shard-dir": args.shard_dir,
            "dataset": args.dataset,
            "checkpoint": args.checkpoint,
            "dry-run-size": args.dry_run_size if args.dry_run_size != 1 else None,
            "n-shards": args.n_shards,
            "shard-index": args.shard_index,
        },
        ("dry-run", "merge-shards", "write-yaml", "retry-failures"),
    )


def write_yaml_results(
    model: Model,
    artifacts: KappaArtifacts,
    merged_run: MergedKappaRun,
) -> None:
    """Evaluate complete predictions and update model YAML metrics/provenance."""
    from matbench_discovery.metrics import phonons as phonon_metrics

    validate_yaml_write_provenance(model, merged_run)
    metrics = phonon_metrics.evaluate_kappa_predictions(artifacts.predictions)
    phonon_metrics.write_metrics_to_yaml(
        model,
        metrics,
        artifacts.pred_file_path,
        run_metadata=merged_run.run_metadata,
        force_file_path=artifacts.force_file_path,
        run_info_path=artifacts.run_info_path,
        replace_pred_file=True,
    )


def validate_yaml_write_provenance(
    model: Model,
    merged_run: MergedKappaRun,
    *,
    canonical_dataset_path: str | None = None,
) -> None:
    """Require official YAML metrics to use the verified protocol and dataset."""
    if merged_run.run_metadata.get("n_failed") != 0:
        raise ValueError("Cannot write kappa YAML metrics from failed records")
    verified_settings = KappaSettings.from_model(model.name)
    if merged_run.manifest.settings != verified_settings:
        raise ValueError(
            "Cannot write kappa YAML metrics from CLI-overridden protocol settings"
        )
    canonical_dataset_path = (
        canonical_dataset_path or DataFiles.phonondb_pbe_103_structures.path
    )
    if merged_run.manifest.dataset_hash != file_sha256(canonical_dataset_path):
        raise ValueError(
            "Cannot write kappa YAML metrics from a noncanonical PhononDB dataset"
        )


def validate_cli_args(
    parser: argparse.ArgumentParser, args: argparse.Namespace
) -> None:
    """Reject incompatible runner modes before loading model dependencies."""
    if args.write_yaml and not args.merge_shards:
        parser.error("--write-yaml is only supported with --merge-shards")
    if args.write_yaml and args.dry_run:
        parser.error("--write-yaml is incompatible with --dry-run")
    if args.merge_shards and args.shard_index is not None:
        parser.error("--shard-index is incompatible with --merge-shards")
    if args.merge_shards and args.retry_failures:
        parser.error("--retry-failures is incompatible with --merge-shards")
    if args.dry_run_size < 1:
        parser.error("--dry-run-size must be positive")
    if not args.dry_run and args.dry_run_size != 1:
        parser.error("--dry-run-size requires --dry-run")


def main(raw_args: Sequence[str] | None = None) -> int:
    """Run, resume, or merge the unified PhononDB kappa benchmark."""
    parser = build_parser()
    args = parser.parse_args(raw_args)
    if args.list_models:
        for model_key in kappa_model_keys():
            spec = CALCULATORS[model_key]
            checkpoint_note = (
                " [requires --checkpoint]" if spec.requires_checkpoint else ""
            )
            dependency_description = (
                f"project={spec.project}"
                if spec.project
                else ", ".join(spec.deps) or "(core deps only)"
            )
            print(f"{model_key}{checkpoint_note}: {dependency_description}")
        return 0
    if not args.model:
        parser.error("--model is required (or pass --list-models)")
    validate_cli_args(parser, args)

    try:
        model_key = resolve_calculator_key(args.model)
        model = Model.from_ref(model_key)
        settings = KappaSettings.from_model(model_key)
    except (TypeError, ValueError) as exc:
        parser.error(f"{exc}, see --list-models")
    if (
        CALCULATORS[model_key].requires_checkpoint
        and not args.checkpoint
        and (not args.merge_shards or args.write_yaml)
    ):
        parser.error(f"{model_key} requires --checkpoint")

    if args.print_cmd:
        command = CALCULATORS[model_key].uv_run_cmd(
            "models/run_kappa.py", *print_cmd_args(args, model_key)
        )
        print(shlex.join(command))
        return 0

    out_dir = args.out_dir or os.path.splitext(model.yaml_path)[0]
    try:
        shard_dir, pred_file_path = resolve_output_paths(
            out_dir=out_dir,
            dry_run=args.dry_run,
            shard_dir=args.shard_dir,
        )
    except ValueError as exc:
        parser.error(str(exc))
    dataset_path = args.dataset or DataFiles.phonondb_pbe_103_structures.path

    if args.merge_shards:
        try:
            atoms_by_id = load_phonondb_atoms(
                dataset_path,
                dry_run=args.dry_run,
                dry_run_size=args.dry_run_size,
            )
            merged_run = merge_kappa_shards(
                shard_dir,
                model_key=model_key,
                expected_material_ids=tuple(atoms_by_id),
                require_103=not args.dry_run,
            )
            if (
                args.n_shards is not None
                and merged_run.manifest.n_shards != args.n_shards
            ):
                raise ValueError(
                    f"Manifest has {merged_run.manifest.n_shards} shards, "
                    f"expected {args.n_shards}"
                )
            if merged_run.manifest.dataset_hash != file_sha256(dataset_path):
                raise ValueError(
                    "Current PhononDB dataset does not match the run manifest"
                )
            if args.write_yaml:
                checkpoint = resolved_artifact_identifier(model, args.checkpoint)
                if merged_run.manifest.checkpoint_hash != checkpoint_digest(checkpoint):
                    raise ValueError(
                        "Current checkpoint does not match the run manifest"
                    )
            artifacts = write_kappa_artifacts(
                merged_run,
                pred_file_path=pred_file_path,
            )
        except (OSError, TypeError, ValueError) as exc:
            parser.error(str(exc))
        print(
            f"Wrote {len(artifacts.predictions):,} predictions "
            f"({artifacts.n_failed:,} failures) to {artifacts.pred_file_path}"
        )
        print(f"Wrote run provenance to {artifacts.run_info_path}")
        if artifacts.force_file_path:
            print(f"Wrote force sets to {artifacts.force_file_path}")
        if args.write_yaml:
            write_yaml_results(model, artifacts, merged_run)
            print(f"Updated kappa metrics and provenance in {model.yaml_path}")
        return int(artifacts.n_failed > 0)

    try:
        n_shards_arg, shard_index_arg, skip_task = effective_shard_args(
            args.n_shards, args.shard_index, dry_run=args.dry_run
        )
        if skip_task:
            print("Skipping redundant dry run outside the first Slurm array task")
            return 0
        n_shards, shard_index = slurm_shard_selection(n_shards_arg, shard_index_arg)
        atoms_by_id = load_phonondb_atoms(
            dataset_path,
            dry_run=args.dry_run,
            dry_run_size=args.dry_run_size,
        )
    except ValueError as exc:
        parser.error(str(exc))

    resolved_device = resolve_device(model_key, args.device)
    calculator = load_calculator(
        model_key,
        device=resolved_device,
        dtype=args.dtype,
        checkpoint=args.checkpoint,
    )
    checkpoint = resolved_artifact_identifier(model, args.checkpoint)
    records = run_kappa_shard(
        calculator=calculator,
        model_key=model_key,
        atoms_by_id=atoms_by_id,
        dataset_path=dataset_path,
        shard_dir=shard_dir,
        shard_index=shard_index,
        n_shards=n_shards,
        settings=settings,
        checkpoint=checkpoint,
        dtype=args.dtype,
        device=resolved_device,
        retry_failures=args.retry_failures,
        dry_run=args.dry_run,
    )
    n_failures = sum(record.error is not None for record in records)
    print(
        f"Completed kappa shard {shard_index + 1}/{n_shards} with "
        f"{len(records):,} records ({n_failures:,} failures) in {shard_dir}"
    )
    return int(n_failures > 0)


if __name__ == "__main__":
    raise SystemExit(main())
