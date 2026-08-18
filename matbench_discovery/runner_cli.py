"""Small helpers shared by dependency-isolated benchmark runners."""

from __future__ import annotations

import glob
import os
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import argparse
    from collections.abc import Mapping, Sequence


def add_common_runner_args(
    parser: argparse.ArgumentParser,
    *,
    dry_run_help: str,
    out_dir_help: str = "Defaults to models/<arch>/<model>",
    list_models_help: str = "Print registered models and exit",
) -> argparse.ArgumentParser:
    """Add the model-selection and dependency-isolation flags every runner shares."""
    parser.add_argument("--model", help="Calculator name, model key, or key alias")
    parser.add_argument("--list-models", action="store_true", help=list_models_help)
    parser.add_argument(
        "--print-cmd",
        action="store_true",
        help="Print the dependency-isolated uv command and exit",
    )
    parser.add_argument("--dry-run", action="store_true", help=dry_run_help)
    parser.add_argument("--out-dir", help=out_dir_help)
    parser.add_argument(
        "--dtype",
        choices=("float64", "float32"),
        default="float64",
        help="Calculator float precision. default float64",
    )
    return parser


def add_shard_args(
    parser: argparse.ArgumentParser, *, merge_shards_help: str, write_yaml_help: str
) -> None:
    """Add the merge/shard flags that sharded runners layer on the common flags."""
    parser.add_argument("--merge-shards", action="store_true", help=merge_shards_help)
    parser.add_argument("--write-yaml", action="store_true", help=write_yaml_help)
    parser.add_argument("--shard-dir", help="Override the resumable shard directory")
    parser.add_argument(
        "--device", choices=("cpu", "cuda"), help="Calculator device. default auto"
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


def resolve_sharded_prefix(
    *,
    default_prefix: str,
    prior_shard_pattern: str,
    task: str,
    shard_dir: str | None,
) -> tuple[str, str]:
    """Select one resumable shard directory and its aligned artifact prefix."""
    selected_shard_dir = shard_dir or f"{default_prefix}-shards"
    if shard_dir is None and not os.path.isdir(selected_shard_dir):
        candidates = sorted(
            path for path in glob.glob(prior_shard_pattern) if os.path.isdir(path)
        )
        if len(candidates) > 1:
            raise ValueError(
                f"Multiple {task} shard directories found: {candidates}. "
                "Select one with --shard-dir."
            )
        if candidates:
            selected_shard_dir = candidates[0]
    selected_shard_dir = os.path.normpath(selected_shard_dir)
    return selected_shard_dir, selected_shard_dir.removesuffix("-shards")


def dependency_run_args(
    args: argparse.Namespace,
    model_key: str,
    value_options: Mapping[str, object],
    flags: Sequence[str],
) -> list[str]:
    """Forward set values and enabled flags to a dependency-isolated command."""
    run_args = ["--model", model_key]
    for option, value in value_options.items():
        if value is not None:
            run_args.extend((f"--{option}", str(value)))
    run_args.extend(
        f"--{flag}" for flag in flags if getattr(args, flag.replace("-", "_"))
    )
    return run_args


def validate_sharded_write_args(
    parser: argparse.ArgumentParser, args: argparse.Namespace
) -> None:
    """Reject write/merge/shard combinations shared by sharded runners."""
    if args.write_yaml and not args.merge_shards:
        parser.error("--write-yaml is only supported with --merge-shards")
    if args.write_yaml and args.dry_run:
        parser.error("--write-yaml is incompatible with --dry-run")
    if args.merge_shards and args.shard_index is not None:
        parser.error("--shard-index is incompatible with --merge-shards")
