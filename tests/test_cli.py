"""Test CLI argument parsing module."""

import pytest

from matbench_discovery import cli
from matbench_discovery.enums import Model, TestSubset
from matbench_discovery.figs import PayloadMode


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
        (
            ["--migrate-model-key", "old-key=new-key"],
            {"migrate_model_key": ("old-key", "new-key")},
            set(),
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


@pytest.mark.parametrize(
    ("selection", "expected"),
    [
        ("models", PayloadMode.targeted),
        ("full_roster", PayloadMode.full_roster),
        ("migrate_provenance", PayloadMode.migrate_provenance),
        ("migrate_model_key", PayloadMode.migrate_model_key),
        (None, None),
    ],
)
def test_payload_mode_is_explicit(
    monkeypatch: pytest.MonkeyPatch,
    selection: str | None,
    expected: PayloadMode | None,
) -> None:
    """Payload mode distinguishes targeted, roster, and migration operations."""
    monkeypatch.setattr(cli, "models_were_explicit", selection == "models")
    monkeypatch.setattr(cli.cli_args, "full_roster", False)
    monkeypatch.setattr(cli.cli_args, "migrate_provenance", False)
    monkeypatch.setattr(cli.cli_args, "migrate_model_key", None)
    if selection not in (None, "models"):
        monkeypatch.setattr(
            cli.cli_args,
            selection,
            ("old-key", "new-key") if selection == "migrate_model_key" else True,
        )
    if expected is None:
        with pytest.raises(SystemExit, match="2"):
            cli.payload_mode()
    else:
        assert cli.payload_mode() == expected
        assert cli.is_full_model_run() is (expected != PayloadMode.targeted)


def test_complete_models_drops_inactive(monkeypatch: pytest.MonkeyPatch) -> None:
    """complete_models drops inactive CLI models even when they appear first."""
    inactive = next(model for model in Model if not model.is_active)
    monkeypatch.setattr(cli.cli_args, "models", [inactive, Model.chgnet_0_3_0])
    assert cli.complete_models() == [Model.chgnet_0_3_0]
