"""Test CLI argument parsing module."""

import pytest

from matbench_discovery import cli
from matbench_discovery.enums import Model, TestSubset


@pytest.mark.parametrize(
    ("args", "expected", "unknown"),
    [
        (
            [],
            {"models": list(Model.active()), "test_subset": TestSubset.uniq_protos},
            set(),
        ),
        (
            ["--models", str(Model.chgnet_0_3_0)],
            {"models": [Model.chgnet_0_3_0]},
            set(),
        ),
        (
            ["--models", "alphanet-v1-mptrj"],
            {"models": [Model.alphanet_v1_mptrj]},
            set(),
        ),
        (
            [
                "--models",
                str(Model.chgnet_0_3_0),
                str(Model.m3gnet),
                "--test-subset",
                "most_stable_10k",
            ],
            {
                "models": [Model.chgnet_0_3_0, Model.m3gnet],
                "test_subset": TestSubset.most_stable_10k,
            },
            set(),
        ),
        (
            [
                "--f=/path/to/kernel.json",
                "--ip=127.0.0.1",
                "--models",
                str(Model.chgnet_0_3_0),
            ],
            {"models": [Model.chgnet_0_3_0]},
            {"--f=/path/to/kernel.json", "--ip=127.0.0.1"},
        ),
    ],
)
def test_cli_parser(
    args: list[str],
    expected: dict[str, TestSubset | list[Model]],
    unknown: set[str],
) -> None:
    """Parse known args; unrecognized Jupyter flags stay in unknown."""
    parsed_args, leftover = cli.cli_parser.parse_known_args(args)
    for key, val in expected.items():
        assert getattr(parsed_args, key) == val
    assert set(leftover) == unknown


@pytest.mark.parametrize(
    ("bad_args", "err_snip"),
    [
        (["--models"], None),
        (["--models", "invalid_model"], "invalid model: invalid_model"),
        (["--test-subset", "invalid_subset"], None),
    ],
)
def test_cli_parser_invalid_args(
    bad_args: list[str], err_snip: str | None, capsys: pytest.CaptureFixture[str]
) -> None:
    """Reject missing/invalid model and subset values."""
    with pytest.raises(SystemExit):
        cli.cli_parser.parse_known_args(bad_args)
    if err_snip:
        error = capsys.readouterr().err
        assert err_snip in error
        assert "None" not in error


def test_shared_payload_test_subset_rejects_model_specific_cohort(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Shared payloads reject the per-model most-stable cohort."""
    monkeypatch.setattr(cli.cli_args, "test_subset", TestSubset.most_stable_10k)
    with pytest.raises(ValueError, match="model-specific"):
        cli.shared_payload_test_subset()


@pytest.mark.parametrize(
    ("models", "is_full"),
    [(list(Model.active()), True), ([Model.chgnet_0_3_0], False)],
)
def test_is_full_model_run(
    monkeypatch: pytest.MonkeyPatch, models: list[Model], is_full: bool
) -> None:
    """is_full_model_run is true only when all active models are selected."""
    monkeypatch.setattr(cli.cli_args, "models", models)
    assert cli.is_full_model_run() is is_full


def test_complete_models_drops_inactive(monkeypatch: pytest.MonkeyPatch) -> None:
    """complete_models drops inactive CLI models even when they appear first."""
    inactive = next(model for model in Model if not model.is_active)
    monkeypatch.setattr(cli.cli_args, "models", [inactive, Model.chgnet_0_3_0])
    assert cli.complete_models() == [Model.chgnet_0_3_0]
