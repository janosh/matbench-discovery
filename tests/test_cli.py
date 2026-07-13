"""Test CLI argument parsing module."""

import pytest

from matbench_discovery.cli import cli_args, cli_parser, shared_payload_test_subset
from matbench_discovery.enums import Model, TestSubset


@pytest.mark.parametrize(
    "args, expected",
    [
        ([], {"models": list(Model.active()), "test_subset": TestSubset.uniq_protos}),
        (["--models", str(Model.chgnet_0_3_0)], {"models": [Model.chgnet_0_3_0]}),
        (
            ["--models", "alphanet-v1-mptrj"],
            {"models": [Model.alphanet_v1_mptrj]},
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
        ),
    ],
)
def test_cli_parser(
    args: list[str], expected: dict[str, TestSubset | list[Model]]
) -> None:
    """Test CLI argument parsing with various inputs."""
    parsed_args, _ = cli_parser.parse_known_args(args)
    for key, val in expected.items():
        assert getattr(parsed_args, key) == val


@pytest.mark.parametrize(
    "bad_args",
    [
        ["--models"],
        ["--models", "invalid_model"],
        ["--test-subset", "invalid_subset"],
    ],
)
def test_cli_parser_invalid_args(
    bad_args: list[str], capsys: pytest.CaptureFixture[str]
) -> None:
    """Test CLI parser raises SystemExit on invalid arguments."""
    with pytest.raises(SystemExit):
        cli_parser.parse_known_args(bad_args)

    if bad_args == ["--models", "invalid_model"]:
        error = capsys.readouterr().err
        assert "invalid model: invalid_model" in error
        assert "None" not in error


def test_cli_parser_jupyter_compat() -> None:
    """Test Jupyter kernel arguments are ignored but preserved in unknown."""
    jupyter_args = [
        "--f=/path/to/kernel.json",
        "--ip=127.0.0.1",
        "--models",
        str(Model.chgnet_0_3_0),
    ]
    args, unknown = cli_parser.parse_known_args(jupyter_args)
    # Our args should be parsed correctly
    assert args.models == [Model.chgnet_0_3_0]
    # Jupyter args should be preserved but ignored
    assert set(unknown) == {"--f=/path/to/kernel.json", "--ip=127.0.0.1"}


def test_shared_payload_test_subset_rejects_model_specific_cohort(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Shared payloads reject the per-model most-stable cohort."""
    monkeypatch.setattr(cli_args, "test_subset", TestSubset.most_stable_10k)
    with pytest.raises(ValueError, match="model-specific"):
        shared_payload_test_subset()


def test_is_full_model_run(monkeypatch: pytest.MonkeyPatch) -> None:
    """is_full_model_run guards multi-model site payloads against filtered runs."""
    from matbench_discovery import cli

    monkeypatch.setattr(cli.cli_args, "models", list(Model.active()))
    assert cli.is_full_model_run() is True

    monkeypatch.setattr(cli.cli_args, "models", [Model.chgnet_0_3_0])
    assert cli.is_full_model_run() is False
