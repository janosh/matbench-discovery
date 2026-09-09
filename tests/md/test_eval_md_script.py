"""Tests for scripts/evals/md.py aggregation."""

import hashlib
import os
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest
import yaml

from matbench_discovery.data import make_file_ref
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
    monkeypatch.setattr(eval_md.cli_args, "md_run_dir", None)
    monkeypatch.setattr(eval_md, "default_md_reference_path", lambda: "ref.h5")
    monkeypatch.setattr(eval_md, "list_reference_systems", lambda _path: systems)
    return eval_md


@pytest.mark.parametrize("use_run_dir", [False, True], ids=["declared", "run-dir"])
@pytest.mark.parametrize(
    ("systems", "column", "problem"),
    [
        (["sysA"], "rdf_error", "missing"),
        (["sysA", "sysB", "sysC", "sysA"], "rdf_error", "duplicate"),
        (["sysA", "sysB", "sysC", "extra"], "rdf_error", "unexpected"),
        (["sysA", "sysB", "sysC"], "unknown_metric", "No recognized MD metric"),
    ],
)
def test_md_evals_preserves_precision_and_rejects_incomplete_coverage(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    systems: list[str],
    column: str,
    problem: str,
    *,
    use_run_dir: bool,
) -> None:
    """Explicit source selection preserves precision and rejects incomplete coverage."""
    model, good_model = Model.mace_mp_0, Model.chgnet_0_3_0
    eval_md = patch_eval_md(tmp_path, monkeypatch, systems=["sysA", "sysB", "sysC"])
    monkeypatch.setattr(eval_md.cli_args, "models", [model, good_model])
    # The generic output flag must never change the selected data source.
    monkeypatch.setattr(eval_md.cli_args, "overwrite", True)

    written_models = []

    def record_write(model: Model, *_args: object, **_kwargs: object) -> None:
        """Record which models passed validation without changing their YAMLs."""
        written_models.append(model)

    monkeypatch.setattr(eval_md.md_metrics, "write_metrics_to_yaml", record_write)

    rdf_error = 2.6901295363762974  # default CSV parsing changes the last binary digit
    invalid_metrics = pd.DataFrame({"system": systems, column: rdf_error})
    # Both sources exist: the unselected source is complete, so choosing it by mistake
    # would pass coverage and reach the YAML writer. Another model succeeds, so a
    # skipped invalid model would incorrectly let the batch report success.
    complete_metrics = pd.DataFrame(
        {"system": ["sysA", "sysB", "sysC"], "rdf_error": [1.0, 2.0, 3.0]}
    )
    combined_csv = tmp_path / "combined.csv.gz"
    (complete_metrics if use_run_dir else invalid_metrics).to_csv(
        combined_csv, index=False
    )
    good_csv = tmp_path / "good.csv.gz"
    complete_metrics.to_csv(good_csv, index=False)
    monkeypatch.setattr(
        Model,
        "md_path",
        property(lambda self: str(combined_csv if self is model else good_csv)),
    )
    md_dir = tmp_path / "models" / os.path.dirname(model.rel_path) / "md-nvt[selected]"
    md_dir.mkdir(parents=True)
    shards = invalid_metrics if use_run_dir else complete_metrics
    for system, frame in shards.groupby("system"):
        frame.to_csv(md_dir / f"{model.name}-md-metrics-{system}.csv.gz", index=False)
    complete_metrics.to_csv(
        md_dir / f"{good_model.name}-md-metrics-all.csv.gz", index=False
    )

    run_dir = str(md_dir) if use_run_dir else None
    monkeypatch.setattr(eval_md.cli_args, "md_run_dir", run_dir)
    resolved = eval_md.resolve_metrics(model, run_dir)
    assert resolved is not None
    assert bool(resolved[1]) is use_run_dir
    assert resolved[0].iloc[0][column] == rdf_error
    assert eval_md.main() == 1
    assert written_models == [good_model]
    assert problem in capsys.readouterr().out
    if resolved[1]:
        assert not (tmp_path / resolved[1]).is_file()


def test_md_recompute_and_import_prediction_file_reference(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Recompute retains source metadata; importing shards writes a fresh artifact."""
    model = Model.mace_mp_0
    eval_md = patch_eval_md(tmp_path, monkeypatch, systems=["sysA"])
    csv_path = tmp_path / "2026-09-12-md-metrics.csv.gz"
    pd.DataFrame({"system": ["sysA"], "rdf_error": [2.5]}).to_csv(csv_path, index=False)
    original_bytes = csv_path.read_bytes()
    ref = make_file_ref(
        str(csv_path),
        url="https://example.org/md.csv.gz",
        size=len(original_bytes),
        md5=hashlib.md5(original_bytes, usedforsecurity=False).hexdigest(),
    )
    metadata = {"model_name": "test model", "metrics": {"md": {"pred_file": ref}}}
    yaml_path = tmp_path / "model.yml"
    yaml_path.write_text(yaml.safe_dump(metadata))
    monkeypatch.setattr(model, "metadata", metadata)
    monkeypatch.setattr(
        type(model), "yaml_path", property(lambda _self: str(yaml_path))
    )
    monkeypatch.setattr(type(model), "md_path", property(lambda _self: str(csv_path)))

    assert eval_md.main() == 0
    result = yaml.safe_load(yaml_path.read_text())["metrics"]["md"]
    assert result == {"pred_file": ref, "rdf_error": 2.5, "n_systems": 1}
    assert csv_path.read_bytes() == original_bytes

    run_dir = tmp_path / "new-run"
    run_dir.mkdir()
    source = pd.DataFrame({"rdf_error": [2.6901295363762974]}, index=["sysA"])
    source.index.name = "system"
    source.to_csv(run_dir / f"{model.name}-md-metrics-sysA.csv.gz")
    monkeypatch.setattr(eval_md.cli_args, "md_run_dir", str(run_dir))
    assert eval_md.main() == 0
    result = yaml.safe_load(yaml_path.read_text())["metrics"]["md"]
    imported = pd.read_csv(
        tmp_path / result["pred_file"]["name"],
        index_col="system",
        float_precision="round_trip",
    )
    pd.testing.assert_frame_equal(imported, source, check_exact=True)
    assert result == {
        "pred_file": {
            "name": f"models/mace/mace-mp-0/{eval_md.today}-md-metrics.csv.gz"
        },
        "rdf_error": 2.6901,
        "n_systems": 1,
    }
    assert csv_path.read_bytes() == original_bytes

    # An empty explicit directory must fail, even with a usable declared artifact.
    monkeypatch.setattr(eval_md.cli_args, "md_run_dir", str(tmp_path / "empty"))
    assert eval_md.main() == 1
    assert yaml.safe_load(yaml_path.read_text())["metrics"]["md"] == result
