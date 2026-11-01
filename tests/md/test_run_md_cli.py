"""Tests for the models/run_md.py command-line interface."""

import shlex
import sys
from pathlib import Path

import pandas as pd
import pytest
from ase.calculators.emt import EMT

from matbench_discovery.calculators import CALCULATORS
from matbench_discovery.enums import Model
from matbench_discovery.metrics import md as md_metrics
from tests.utils import import_repo_script

RUN_MD = import_repo_script("run_md", "models/run_md.py")
DEFAULT_MD_OPTIONS = shlex.split(
    "--n-steps 80000 --time-step-fs 0.25 --record-interval 10 --seed 0 --dtype float64"
)
CUSTOM_MD_OPTIONS = shlex.split(
    "--out-dir /tmp/md-out --ref-file /tmp/ref.h5 "
    "--private-ref-file /tmp/private-ref.h5 --n-steps 123 --time-step-fs 0.5 "
    "--record-interval 4 --seed 7 --dtype float32 --dry-run --systems bulkAu bulkCu"
)
WRITE_MD_OPTIONS = [*DEFAULT_MD_OPTIONS, "--write-yaml"]
REF_SYSTEMS = ("bulkAu", "bulkCu", "bulkSi")


def fake_load_calculator(
    *_args: object, **_kwargs: object
) -> tuple[EMT, dict[str, object]]:
    """Return a lightweight calculator and its runtime for CLI tests."""
    return EMT(), {"hardware": "CPU (test)"}


def test_run_md_cli_write_yaml_skips_non_submission_model(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`--write-yaml` for debug model emt must skip YAML write gracefully."""
    # avoid the real calculator + rollout pipeline; return metrics calc_md_metrics reads
    captured_kwargs: dict[str, object] = {}

    def fake_run_md_benchmark(**kwargs: object) -> pd.DataFrame:
        """Capture runner kwargs while returning minimal aggregateable metrics."""
        captured_kwargs.update(kwargs)
        return pd.DataFrame({"rdf_error": [1.0], "vdos_error": [2.0]})

    monkeypatch.setattr(RUN_MD, "load_calculator_with_runtime", fake_load_calculator)
    monkeypatch.setattr(RUN_MD, "run_md_benchmark", fake_run_md_benchmark)
    monkeypatch.setenv("MBD_MD_PRIVATE_REF_FILE", "/private/ref.h5")
    monkeypatch.setattr(sys, "argv", ["run_md", "--model", "emt", "--write-yaml"])

    assert RUN_MD.main() == 0  # would raise ValueError without the non-enum guard
    assert captured_kwargs["private_ref_file"] == "/private/ref.h5"
    assert captured_kwargs["runtime"] == {"hardware": "CPU (test)"}


@pytest.mark.parametrize(
    ("systems", "problem"),
    [
        (REF_SYSTEMS[:-1], "1 missing e.g. ['bulkSi']"),
        ((*REF_SYSTEMS, "bulkCu"), "duplicate systems ['bulkCu']"),
    ],
    ids=["missing", "duplicated"],
)
def test_run_md_cli_write_yaml_rejects_incomplete_coverage(
    systems: tuple[str, ...],
    problem: str,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """`--write-yaml` refuses to publish metrics.md without exact system coverage."""
    yaml_path = Model["mace_mp_0"].yaml_path
    with open(yaml_path, encoding="utf-8") as file:
        yaml_before = file.read()
    yaml_writes: list[object] = []
    monkeypatch.setattr(RUN_MD, "load_calculator_with_runtime", fake_load_calculator)
    monkeypatch.setattr(
        RUN_MD,
        "run_md_benchmark",
        lambda **_kwargs: pd.DataFrame(
            {"rdf_error": 1.0, "vdos_error": 2.0}, index=list(systems)
        ),
    )
    monkeypatch.setattr(RUN_MD, "list_reference_systems", lambda _path: REF_SYSTEMS)
    monkeypatch.setattr(RUN_MD, "ROOT", str(tmp_path))  # catch any CSV write
    monkeypatch.setattr(
        md_metrics,
        "write_metrics_to_yaml",
        lambda *args, **_kw: yaml_writes.append(args),
    )
    monkeypatch.setattr(sys, "argv", ["run_md", "--model", "mace_mp_0", "--write-yaml"])

    with pytest.raises(
        ValueError, match=r"Not writing metrics\.md for mace_mp_0"
    ) as exc:
        RUN_MD.main()
    assert problem in str(exc.value)
    assert yaml_writes == []
    assert list(tmp_path.iterdir()) == []
    with open(yaml_path, encoding="utf-8") as file:
        assert file.read() == yaml_before


@pytest.mark.parametrize(
    ("argv", "expected_cmd"),
    [
        (
            ["run_md", "--model", "mace_mp_0", "--write-yaml", "--systems", "bulkAu"],
            None,
        ),
        (
            ["run_md", "--print-cmd", "--model", "mace-mp-0", *CUSTOM_MD_OPTIONS],
            CALCULATORS["mace_mp_0"].uv_run_cmd(
                "models/run_md.py", "--model", "mace_mp_0", *CUSTOM_MD_OPTIONS
            ),
        ),
        (
            ["run_md", "--print-cmd", "--model", "mace-mp-0", *WRITE_MD_OPTIONS],
            CALCULATORS["mace_mp_0"].uv_run_cmd(
                "models/run_md.py", "--model", "mace_mp_0", *WRITE_MD_OPTIONS
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
        assert "--write-yaml needs a full run" in capsys.readouterr().err
    else:
        assert RUN_MD.main() == 0
        assert shlex.split(capsys.readouterr().out) == expected_cmd
