"""Tests for the models/run_md.py command-line interface."""

import shlex
import sys

import pandas as pd
import pytest
from ase.calculators.emt import EMT

from matbench_discovery.calculators import CALCULATORS
from tests.utils import import_repo_script

RUN_MD = import_repo_script("run_md", "models/run_md.py")


def test_run_md_cli_write_yaml_skips_non_submission_model(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`--write-yaml` for debug model emt must skip YAML write gracefully."""
    # avoid the real calculator + rollout pipeline; return metrics calc_md_metrics reads
    captured_kwargs: dict[str, object] = {}

    def fake_load_calculator(*_args: object, **_kwargs: object) -> EMT:
        """Return a lightweight calculator for the CLI test."""
        return EMT()

    def fake_run_md_benchmark(**kwargs: object) -> pd.DataFrame:
        """Capture runner kwargs while returning minimal aggregateable metrics."""
        captured_kwargs.update(kwargs)
        return pd.DataFrame({"rdf_error": [1.0], "vdos_error": [2.0]})

    monkeypatch.setattr(RUN_MD, "load_calculator", fake_load_calculator)
    monkeypatch.setattr(RUN_MD, "run_md_benchmark", fake_run_md_benchmark)
    monkeypatch.setenv("MBD_MD_PRIVATE_REF_FILE", "/private/ref.h5")
    monkeypatch.setattr(sys, "argv", ["run_md", "--model", "emt", "--write-yaml"])

    assert RUN_MD.main() == 0  # would raise ValueError without the non-enum guard
    assert captured_kwargs["private_ref_file"] == "/private/ref.h5"


@pytest.mark.parametrize(
    ("argv", "expected_cmd"),
    [
        (
            ["run_md", "--model", "mace_mp_0", "--write-yaml", "--systems", "bulkAu"],
            None,
        ),
        (
            ["run_md", "--print-cmd", "--model", "mace-mp-0"],
            CALCULATORS["mace_mp_0"].uv_run_cmd(
                "models/run_md.py", "--model", "mace_mp_0"
            ),
        ),
    ],
)
def test_run_md_cli_validation(
    argv: list[str],
    expected_cmd: list[str] | None,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """CLI rejects partial YAML writes and prints full canonical uv commands."""
    monkeypatch.setattr(sys, "argv", argv)

    if expected_cmd is None:
        with pytest.raises(SystemExit):  # parser.error exits before rollout pipeline
            RUN_MD.main()
    else:
        assert RUN_MD.main() == 0
        assert shlex.split(capsys.readouterr().out) == expected_cmd
