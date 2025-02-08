"""Test CLI argument parsing module."""

import sys
from unittest.mock import patch

import pytest
from pymatviz.enums import Key

from matbench_discovery.cli import cli_args, cli_parser
from matbench_discovery.enums import Model, TestSubset


@pytest.mark.parametrize(
    "args, expected",
    [
        ([], {"models": list(Model), "test_subset": TestSubset.uniq_protos}),
        (["--models", "chgnet"], {"models": [Model.chgnet]}),
        (
            ["--models", "chgnet", "m3gnet", "--test-subset", "most_stable_10k"],
            {
                "models": [Model.chgnet, Model.m3gnet],
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
    with patch.object(sys, "argv", ["script.py", *args]):
        parsed_args, _ = cli_parser.parse_known_args()
        for key, val in expected.items():
            assert getattr(parsed_args, key) == val


@pytest.mark.parametrize(
    "bad_args",
    [
        ["--models", "invalid_model"],
        ["--test-subset", "invalid_subset"],
        ["--energy-type", "invalid_type"],
    ],
)
def test_cli_parser_invalid_args(bad_args: list[str]) -> None:
    """Test CLI parser raises SystemExit on invalid arguments."""
    with pytest.raises(SystemExit), patch.object(sys, "argv", ["script.py", *bad_args]):
        cli_parser.parse_known_args()


def test_cli_parser_jupyter_compat() -> None:
    """Test Jupyter kernel arguments are ignored but preserved in unknown."""
    jupyter_args = ["--f=/path/to/kernel.json", "--ip=127.0.0.1", "--models", "chgnet"]
    with patch.object(sys, "argv", ["script.py", *jupyter_args]):
        args, unknown = cli_parser.parse_known_args()
        # Our args should be parsed correctly
        assert args.models == [Model.chgnet]
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
