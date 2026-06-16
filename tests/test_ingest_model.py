"""Tests for the submission-ingestion checklist (scripts/ingest_model.py).

The checklist reads model YAMLs through the Model enum metadata API (replacing the
old grep/sed bash recipe), so it can be tested against the real repo metadata
without running any evals.
"""

import pytest

import scripts.ingest_model as ingest
from matbench_discovery.enums import Model


def msgs(checks: "ingest.Checklist", status: str) -> list[str]:
    return [msg for stat, msg in checks.results if stat == status]


def test_complete_force_model_passes_checklist() -> None:
    """A complete UIP submission (mace-mpa-0) passes every applicable check."""
    checks = ingest.Checklist()
    energy_only = ingest.check_submission(Model.mace_mpa_0, checks)
    assert energy_only is False
    assert msgs(checks, ingest.FAIL) == []
    # force-based checks must actually run, not skip
    passed = msgs(checks, ingest.PASS)
    assert any("geo_opt" in msg for msg in passed)
    assert any("Phonon" in msg for msg in passed)


def test_energy_only_model_skips_force_tasks() -> None:
    """targets=E models skip geo-opt/phonons/diatomics instead of failing."""
    e_only = next((mdl for mdl in Model if mdl.metadata.get("targets") == "E"), None)
    if e_only is None:
        pytest.skip("no energy-only model in registry")
    checks = ingest.Checklist()
    assert ingest.check_submission(e_only, checks) is True
    assert msgs(checks, ingest.FAIL) == []
    skips = msgs(checks, ingest.SKIP)
    assert sum("skipped (targets=E" in msg for msg in skips) >= 4


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
            msg for msg in msgs(checks, ingest.FAIL) if "test script" not in msg
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


def test_run_payload_refresh_models_flag(monkeypatch: pytest.MonkeyPatch) -> None:
    """Payload refresh passes --models for single-model runs (merge mode) and omits
    it for full regens; the final command runs the payload shape tests either way.
    """
    calls: list[tuple[str, ...]] = []

    def fake_run_cmd(*cmd: str) -> bool:
        calls.append(cmd)
        return True

    monkeypatch.setattr(ingest, "run_cmd", fake_run_cmd)
    checks = ingest.Checklist()
    ingest.run_payload_refresh(checks, model=Model.mace_mpa_0)
    *script_calls, test_call = calls
    assert len(script_calls) == len(ingest.PAYLOAD_SCRIPTS)
    assert any(
        "single_model_per_element_errors.py" in " ".join(cmd) for cmd in script_calls
    )
    for cmd in script_calls:
        assert cmd[-2:] == ("--models", "mace_mpa_0")
    assert "pytest" in test_call
    assert checks.n_failed == 0

    calls.clear()
    ingest.run_payload_refresh(ingest.Checklist())
    assert all("--models" not in cmd for cmd in calls)


def test_map_yaml_paths() -> None:
    """--map-yaml-paths maps PR-diff YAML paths to enum names, fails on unknowns."""
    path = Model.mace_mpa_0.yaml_path.split("/models/")[-1]
    assert ingest.map_yaml_paths([f"models/{path}", f"models/{path}"]) == [
        "mace_mpa_0"  # deduped
    ]
    with pytest.raises(SystemExit, match="No Model enum entry"):
        ingest.map_yaml_paths(["models/new-arch/unregistered.yml"])
