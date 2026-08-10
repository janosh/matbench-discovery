"""Tests for the per-model sweep shared by scripts/evals."""

import pytest

from matbench_discovery.enums import Model
from scripts.evals import evaluate_models

MODEL_A, MODEL_B = Model.mace_mp_0, Model.chgnet_0_3_0
Outcome = str | Exception | None  # skip reason, raised error, or success


@pytest.mark.parametrize(
    ("models", "outcomes", "expected_code"),
    [
        ([], {}, 0),  # selecting nothing is not a failure
        ([MODEL_A], {MODEL_A: None}, 0),
        ([MODEL_A], {MODEL_A: "no artifact"}, 1),
        ([MODEL_A], {MODEL_A: ValueError("boom")}, 1),
        ([MODEL_A, MODEL_B], {MODEL_A: None, MODEL_B: "no artifact"}, 0),
        ([MODEL_A, MODEL_B], {MODEL_A: None, MODEL_B: OSError("boom")}, 1),
    ],
    ids=["empty", "ok", "all-skipped", "all-failed", "partial-skip", "partial-fail"],
)
def test_evaluate_models_exit_code(
    models: list[Model], outcomes: dict[Model, Outcome], expected_code: int
) -> None:
    """Errors fail the sweep even alongside successes; skips and empty runs don't."""

    def evaluate_one(model: Model) -> str | None:
        outcome = outcomes[model]
        if isinstance(outcome, Exception):
            raise outcome
        return outcome

    assert evaluate_models("test", models, evaluate_one) == expected_code


def test_evaluate_models_counts_errors_apart_from_skips(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The summary must not report an errored model as skipped."""

    def evaluate_one(model: Model) -> str | None:
        if model is MODEL_B:
            raise KeyError("boom")
        return "no artifact"

    assert evaluate_models("test", [MODEL_A, MODEL_B], evaluate_one) == 1
    assert "1 skipped, 1 failed" in capsys.readouterr().out
