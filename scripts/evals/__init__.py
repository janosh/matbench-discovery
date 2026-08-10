"""Evaluation scripts for Matbench Discovery."""

import traceback
from collections.abc import Callable, Sequence

from matbench_discovery.enums import Model


def evaluate_models(
    task: str, models: Sequence[Model], evaluate_one: Callable[[Model], str | None]
) -> int:
    """Run one evaluator per model and report success."""
    print(f"Evaluating {task} metrics for {len(models)} model(s)...")
    n_success = 0
    for model in models:
        try:
            skip_reason = evaluate_one(model)
        except (ValueError, OSError, KeyError) as exc:
            print(f"\tError processing {model.label}: {exc}")
            traceback.print_exc()  # unexpected, so keep the full context
            continue
        if skip_reason is None:
            n_success += 1
        else:
            print(f"Skipping {model.label}: {skip_reason}")

    n_skipped = len(models) - n_success
    if n_success == 0:
        print(f"\nNo models evaluated successfully ({n_skipped} skipped)")
        return 1
    print(f"\nSuccessfully evaluated {n_success} model(s), {n_skipped} skipped")
    return 0
