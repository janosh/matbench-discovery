"""Test CLI argument parsing module."""

import pytest
from pymatviz.enums import Key

from matbench_discovery.cli import cli_args, cli_parser
from matbench_discovery.enums import Model, TestSubset


@pytest.mark.parametrize(
    "args, expected",
    [
        ([], {"models": list(Model.active()), "test_subset": TestSubset.uniq_protos}),
        (["--models", str(Model.chgnet_030)], {"models": [Model.chgnet_030]}),
        (
            ["--models", "alphanet-mptrj"],
            {"models": [Model.alphanet_mptrj]},
        ),
        (
            [
                "--models",
                str(Model.chgnet_030),
                str(Model.m3gnet_ms),
                "--test-subset",
                "most_stable_10k",
            ],
            {
                "models": [Model.chgnet_030, Model.m3gnet_ms],
                "test_subset": TestSubset.most_stable_10k,
            },
        ),
        (
            ["--energy-type", "e_form", "--show-non-compliant"],
            {"energy_type": Key.e_form, "show_non_compliant": True},
        ),
    ],
)
def test_cli_parser(
    args: list[str], expected: dict[str, str | bool | TestSubset | list[Model]]
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
        ["--energy-type", "invalid_type"],
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
        str(Model.chgnet_030),
    ]
    args, unknown = cli_parser.parse_known_args(jupyter_args)
    # Our args should be parsed correctly
    assert args.models == [Model.chgnet_030]
    # Jupyter args should be preserved but ignored
    assert set(unknown) == {"--f=/path/to/kernel.json", "--ip=127.0.0.1"}


def test_cli_args_global() -> None:
    """Test global cli_args is properly initialized."""
    assert cli_args is not None
    # Test all expected attributes are present and of correct type
    assert isinstance(cli_args.models, list)
    assert isinstance(cli_args.test_subset, TestSubset)
    assert isinstance(cli_args.energy_type, str)
    assert isinstance(cli_args.show_non_compliant, bool)
    assert isinstance(cli_args.use_full_rows, bool)
    assert isinstance(cli_args.update_existing, bool)


def test_browser_renderers_never_steal_focus() -> None:
    """Figures opened as browser tabs must not autoraise (switch screen focus)."""
    import plotly.io as pio

    for renderer_name in ("browser", "chrome", "chromium", "firefox"):
        assert pio.renderers[renderer_name].autoraise is False, renderer_name


def test_is_full_model_run(monkeypatch: pytest.MonkeyPatch) -> None:
    """is_full_model_run guards multi-model site payloads against filtered runs."""
    from matbench_discovery import cli

    monkeypatch.setattr(cli.cli_args, "models", list(Model.active()))
    assert cli.is_full_model_run() is True

    monkeypatch.setattr(cli.cli_args, "models", [Model.chgnet_030])
    assert cli.is_full_model_run() is False
