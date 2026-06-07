"""Tests for the submission-ingestion checklist (scripts/ingest_model.py).

The checklist reads model YAMLs through the Model enum metadata API (replacing the
old grep/sed bash recipe), so it can be tested against the real repo metadata
without running any evals.
"""

import importlib.util
import sys
from pathlib import Path

import pytest

from matbench_discovery.enums import Model


def import_script(name: str):  # noqa: ANN201
    """Import a module from scripts/ (not a package; can't live in conftest since
    `from tests.conftest import ...` resolves to pymatviz's tests package).
    """
    path = Path(__file__).parent.parent / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None  # narrow ModuleSpec | None for type checker
    assert spec.loader is not None  # narrow Loader | None for type checker
    sys.modules[name] = module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ingest = import_script("ingest_model")


def results_by_status(checks: "ingest.Checklist") -> dict[str, list[str]]:
    out: dict[str, list[str]] = {ingest.PASS: [], ingest.FAIL: [], ingest.SKIP: []}
    for status, msg in checks.results:
        out[status].append(msg)
    return out


def test_complete_force_model_passes_checklist() -> None:
    """A complete UIP submission (mace-mpa-0) passes every applicable check."""
    checks = ingest.Checklist()
    energy_only = ingest.check_submission(Model.mace_mpa_0, checks)
    assert energy_only is False
    by_status = results_by_status(checks)
    assert by_status[ingest.FAIL] == []
    # force-based checks must actually run, not skip
    assert any("geo_opt" in msg for msg in by_status[ingest.PASS])
    assert any("Phonon" in msg for msg in by_status[ingest.PASS])


def test_energy_only_model_skips_force_tasks() -> None:
    """targets=E models skip geo-opt/phonons/diatomics instead of failing."""
    e_only = next((mdl for mdl in Model if mdl.metadata.get("targets") == "E"), None)
    if e_only is None:
        pytest.skip("no energy-only model in registry")
    checks = ingest.Checklist()
    assert ingest.check_submission(e_only, checks) is True
    by_status = results_by_status(checks)
    assert by_status[ingest.FAIL] == [], by_status[ingest.FAIL]
    assert sum("skipped (targets=E" in msg for msg in by_status[ingest.SKIP]) >= 4


def test_all_active_models_have_required_metadata() -> None:
    """Every active model has the YAML + discovery pred URL that evals and the site
    depend on. Test-script checks are excluded: pre-kappa-era models are
    grandfathered without test_*_kappa.py and the checklist only gates NEW
    submissions on those.
    """
    failures: dict[str, list[str]] = {}
    for model in Model.active():
        checks = ingest.Checklist()
        ingest.check_submission(model, checks)
        hard_fails = [
            msg
            for status, msg in checks.results
            if status == ingest.FAIL and "test script" not in msg
        ]
        if hard_fails:
            failures[model.name] = hard_fails
    assert not failures, failures


def test_checklist_summary_counts() -> None:
    checks = ingest.Checklist()
    checks.ok("a")
    checks.fail("b")
    checks.skip("c")
    checks.skip("d")
    assert checks.n_failed == 1
    assert "Passed:  1" in checks.summary()
    assert "Failed:  1" in checks.summary()
    assert "Skipped: 2" in checks.summary()


def test_cli_rejects_unknown_model_and_missing_args() -> None:
    with pytest.raises(SystemExit, match="2"):
        ingest.main(["definitely-not-a-model"])
    with pytest.raises(SystemExit, match="2"):
        ingest.main([])  # model required unless --payloads-only


def test_cli_archive_requires_figshare_token(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv("FIGSHARE_TOKEN", raising=False)
    with pytest.raises(SystemExit, match="2"):
        ingest.main(["mace-mpa-0", "--archive"])


def test_map_yaml_paths() -> None:
    """--map-yaml-paths maps PR-diff YAML paths to enum names, fails on unknowns."""
    path = Model.mace_mpa_0.yaml_path.split("/models/")[-1]
    assert ingest.map_yaml_paths([f"models/{path}", f"models/{path}"]) == [
        "mace_mpa_0"  # deduped
    ]
    with pytest.raises(SystemExit, match="No Model enum entry"):
        ingest.map_yaml_paths(["models/new-arch/unregistered.yml"])
