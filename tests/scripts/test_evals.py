"""Tests for the per-model sweep shared by scripts/evals."""

import pytest

from matbench_discovery.enums import Model
from scripts.evals import evaluate_models

MODEL_A, MODEL_B = Model.mace_mp_0, Model.chgnet_0_3_0
Outcome = str | Exception | None  # skip reason, raised error, or success


@pytest.mark.parametrize(
    ("outcomes", "expected_code", "expected_tally"),
    [
        ({}, 0, "0 skipped, 0 failed"),  # selecting nothing is not a failure
        ({MODEL_A: None}, 0, "0 skipped, 0 failed"),
        ({MODEL_A: "no artifact"}, 1, "1 skipped, 0 failed"),
        ({MODEL_A: ValueError("boom")}, 1, "0 skipped, 1 failed"),
        ({MODEL_A: None, MODEL_B: "no artifact"}, 0, "1 skipped, 0 failed"),
        ({MODEL_A: None, MODEL_B: OSError("boom")}, 1, "0 skipped, 1 failed"),
        ({MODEL_A: "no artifact", MODEL_B: KeyError("boom")}, 1, "1 skipped, 1 failed"),
    ],
    ids=["empty", "ok", "skip", "error", "ok-skip", "ok-error", "skip-error"],
)
def test_evaluate_models_exit_code(
    outcomes: dict[Model, Outcome],
    expected_code: int,
    expected_tally: str,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Exit codes and summaries distinguish successful, skipped, and failed models."""

    def evaluate_one(model: Model) -> str | None:
        outcome = outcomes[model]
        if isinstance(outcome, Exception):
            raise outcome
        return outcome

    assert evaluate_models("test", list(outcomes), evaluate_one) == expected_code
    assert expected_tally in capsys.readouterr().out
