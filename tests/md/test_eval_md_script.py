"""Tests for scripts/evals/md.py aggregation."""

import os
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

from matbench_discovery.enums import Model
from tests.utils import import_repo_script


def patch_eval_md(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    systems: list[str],
    model: Model = Model.mace_mp_0,
) -> ModuleType:
    """Import scripts/evals/md.py and patch it to use a tiny test reference set."""
    eval_md = import_repo_script("eval_md", "scripts/evals/md.py")
    monkeypatch.setattr(eval_md, "ROOT", str(tmp_path))
    monkeypatch.setattr(eval_md.cli_args, "models", [model])
    monkeypatch.setattr(eval_md, "default_md_reference_path", lambda: "ref.h5")
    monkeypatch.setattr(eval_md, "list_reference_systems", lambda _path: systems)
    return eval_md


@pytest.mark.parametrize("fallback", [False, True], ids=["per-system", "fallback"])
def test_md_evals_skips_incomplete_coverage(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *, fallback: bool
) -> None:
    """Per-system and fallback CSVs must both reject incomplete coverage."""
    model = Model.mace_mp_0
    eval_md = patch_eval_md(tmp_path, monkeypatch, systems=["sysA", "sysB", "sysC"])

    def _fail_write(*_args: object, **_kwargs: object) -> None:
        """Fail if incomplete coverage reaches the YAML writer."""
        raise AssertionError("must not write YAML on incomplete coverage")

    monkeypatch.setattr(eval_md.md_metrics, "write_metrics_to_yaml", _fail_write)

    incomplete_metrics = pd.DataFrame({"system": ["sysA"], "rdf_error": [1.0]})
    if fallback:
        # md_path fallback: submitted combined CSV with no per-system files
        incomplete_csv = tmp_path / "combined.csv.gz"  # only 1 of 3 reference systems
        incomplete_metrics.to_csv(incomplete_csv, index=False)
        monkeypatch.setattr(
            type(model), "md_path", property(lambda _self: str(incomplete_csv))
        )
    else:
        arch_dir = os.path.dirname(model.rel_path)
        md_dir = tmp_path / "models" / arch_dir / "2026-06-14-md-nvt"
        md_dir.mkdir(parents=True)
        incomplete_metrics.to_csv(
            md_dir / f"{model.name}-md-metrics-sysA.csv.gz", index=False
        )

    assert eval_md.main() == 1  # missing sysB/sysC -> skip, exit 1, no YAML written


def test_md_evals_handles_placeholder_md_metrics(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A model with complete per-system coverage but a 'not available' placeholder for
    metrics.md (a string, not a dict) must aggregate without crashing on md_yaml.get().
    """
    model = Model.mace_mp_0
    monkeypatch.setattr(
        type(model), "metrics", property(lambda _self: {"md": "not available"})
    )
    eval_md = patch_eval_md(tmp_path, monkeypatch, systems=["sysA"])

    writes: list[dict[str, object]] = []
    monkeypatch.setattr(
        eval_md.md_metrics,
        "write_metrics_to_yaml",
        lambda *_args, **kwargs: writes.append(kwargs),
    )

    arch_dir = os.path.dirname(model.rel_path)
    md_dir = tmp_path / "models" / arch_dir / "2026-06-14-md-nvt"
    md_dir.mkdir(parents=True)
    pd.DataFrame({"system": ["sysA"], "rdf_error": [1.0], "vdos_error": [2.0]}).to_csv(
        md_dir / f"{model.name}-md-metrics-sysA.csv.gz", index=False
    )

    assert eval_md.main() == 0  # placeholder md must not raise AttributeError
    assert len(writes) == 1
    assert writes[0]["pred_file_url"] is None  # no url from a placeholder md section
