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
                "full_test_set",
            ],
            {
                "models": [Model.chgnet_0_3_0, Model.m3gnet],
                "test_subset": TestSubset.full_test_set,
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
    expected: dict[str, object],
    unknown: set[str],
) -> None:
    """Parse known args; unrecognized Jupyter flags stay in unknown."""
    parsed_args, leftover = cli.cli_parser.parse_known_args(args)
    for key, val in expected.items():
        assert getattr(parsed_args, key) == val
    assert set(leftover) == unknown
    explicit_models = any(arg.partition("=")[0] == "--models" for arg in args)
    assert (parsed_args.models is not cli.models_arg.default) is explicit_models


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


@pytest.mark.parametrize("models_explicit", [False, True])
@pytest.mark.parametrize("full_roster", [False, True])
def test_payload_scope_is_explicit(
    monkeypatch: pytest.MonkeyPatch,
    models_explicit: bool,
    full_roster: bool,
) -> None:
    """Payload generation requires exactly one targeted or full-roster scope."""
    monkeypatch.setattr(cli, "models_were_explicit", models_explicit)
    monkeypatch.setattr(cli.cli_args, "full_roster", full_roster)
    if models_explicit == full_roster:
        with pytest.raises(SystemExit, match="2"):
            cli.is_full_model_run()
    else:
        assert cli.is_full_model_run() is full_roster


def test_complete_models_drops_inactive(monkeypatch: pytest.MonkeyPatch) -> None:
    """complete_models drops inactive CLI models even when they appear first."""
    inactive = next(model for model in Model if not model.is_active)
    monkeypatch.setattr(cli.cli_args, "models", [inactive, Model.chgnet_0_3_0])
    assert cli.complete_models() == [Model.chgnet_0_3_0]
