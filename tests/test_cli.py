"""Test CLI argument parsing module."""

import sys
from types import ModuleType

import pytest

from matbench_discovery import cli
from matbench_discovery.enums import Model, TestSubset
from tests.utils import patch_model_metadata


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
            ["--models", "alphanet-v1-mptrj", "--md-run-dir", "models/run-2026-09-12"],
            {
                "models": [Model.alphanet_v1_mptrj],
                "md_run_dir": "models/run-2026-09-12",
            },
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
        (["--md-run-dir"], None),
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
    models = [Model.chgnet_0_3_0] if models_explicit else cli.models_arg.default
    monkeypatch.setattr(cli.cli_args, "models", models)
    monkeypatch.setattr(cli.cli_args, "full_roster", full_roster)
    if models_explicit == full_roster:
        with pytest.raises(SystemExit, match="2"):
            cli.is_full_model_run()
    else:
        assert cli.is_full_model_run() is full_roster


def test_complete_models(monkeypatch: pytest.MonkeyPatch) -> None:
    """complete_models keeps only active CLI models with a discovery pred_file."""
    pred_metrics = {"discovery": {"pred_file": {"name": "models/fake/preds.csv.gz"}}}
    metadata_by_model: dict[Model, dict[str, object]] = {
        # inactive listed first must not shadow later complete models
        Model.alignn: {"lifecycle": "superseded", "metrics": pred_metrics},
        Model.chgnet_0_3_0: {"lifecycle": "active", "metrics": pred_metrics},
        Model.m3gnet: {"lifecycle": "active", "metrics": None},  # bare `metrics:`
        Model.mace_mp_0: {"lifecycle": "active", "metrics": {"discovery": None}},
    }
    for model, metadata in metadata_by_model.items():
        patch_model_metadata(monkeypatch, model, metadata)
    monkeypatch.setattr(cli.cli_args, "models", list(metadata_by_model))
    assert cli.complete_models() == [Model.chgnet_0_3_0]


@pytest.mark.parametrize(
    ("argv", "jupyter", "expected_dry_run"),
    [
        (["-n"], False, True),
        (["--dryrun"], False, None),  # typo exits instead of silently defaulting
        (["--f=/path/to/kernel.json"], True, False),  # Jupyter kernel flags pass
    ],
)
def test_parse_cli_args_is_strict(
    monkeypatch: pytest.MonkeyPatch,
    argv: list[str],
    jupyter: bool,
    expected_dry_run: bool | None,
) -> None:
    """Shared CLI parsing rejects unknown flags outside Jupyter and caches args."""
    monkeypatch.setattr(cli, "cli_args", cli.cli_args)  # restore cache afterwards
    if jupyter:
        monkeypatch.setitem(sys.modules, "ipykernel", ModuleType("ipykernel"))
    else:
        monkeypatch.delitem(sys.modules, "ipykernel", raising=False)
    if expected_dry_run is None:
        with pytest.raises(SystemExit, match="2"):
            cli.parse_cli_args(argv)
        return
    args = cli.parse_cli_args(argv)
    assert args.dry_run is expected_dry_run
    assert cli.cli_args is args
