"""Small helpers shared by dependency-isolated benchmark runners."""

from __future__ import annotations

import glob
import os
import shlex
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import argparse
    from collections.abc import Mapping, Sequence


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


def print_dependency_command(command: Sequence[str]) -> None:
    """Print a dependency-isolated command with shell-safe quoting."""
    print(shlex.join(command))


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
