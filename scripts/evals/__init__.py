"""Evaluation scripts for Matbench Discovery."""

import traceback
from collections.abc import Callable, Sequence

from matbench_discovery.enums import Model


def evaluate_models(
    task: str, models: Sequence[Model], evaluate_one: Callable[[Model], str | None]
) -> int:
    """Run one evaluator per model, returning 0 only if nothing failed.

    Evaluators return a reason string to skip an expected gap. Selecting no models is
    not a failure, but an unexpected error is, even when other models succeed.
    """
    print(f"Evaluating {task} metrics for {len(models)} model(s)...")
    n_success = n_failed = 0
    for model in models:
        try:
            skip_reason = evaluate_one(model)
        except (ValueError, OSError, KeyError) as exc:
            print(f"\tError processing {model.label}: {exc}")
            traceback.print_exc()  # unexpected, so keep the full context
            n_failed += 1
            continue
        if skip_reason is None:
            n_success += 1
        else:
            print(f"Skipping {model.label}: {skip_reason}")

    tally = f"{len(models) - n_success - n_failed} skipped, {n_failed} failed"
    if n_success == 0:
        print(f"\nNo models evaluated successfully ({tally})")
    else:
        print(f"\nSuccessfully evaluated {n_success} model(s), {tally}")
    return 1 if n_failed or (len(models) > 0 and n_success == 0) else 0
