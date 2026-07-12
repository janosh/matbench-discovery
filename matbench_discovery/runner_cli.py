"""Small helpers shared by dependency-isolated benchmark runners."""

from __future__ import annotations

import glob
import os
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
