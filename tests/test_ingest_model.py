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


@pytest.fixture
def run_cmd_calls(monkeypatch: pytest.MonkeyPatch) -> list[tuple[str, ...]]:
    """Patch ingest.run_cmd and return captured commands."""
    calls: list[tuple[str, ...]] = []

    def fake_run_cmd(*cmd: str) -> bool:
        calls.append(cmd)
        return True

    monkeypatch.setattr(ingest, "run_cmd", fake_run_cmd)
    return calls


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
    assert e_only is not None
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


@pytest.mark.parametrize("argv", [["definitely-not-a-model"], []])
def test_cli_rejects_unknown_model_and_missing_args(argv: list[str]) -> None:
    """Unknown or missing model args fail argparse validation."""
    with pytest.raises(SystemExit, match="2"):
        ingest.main(argv)


def test_cli_archive_requires_figshare_token(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("FIGSHARE_TOKEN", raising=False)
    with pytest.raises(SystemExit, match="2"):
        ingest.main(["mace-mpa-0", "--archive"])


def test_run_payload_refresh_models_flag(run_cmd_calls: list[tuple[str, ...]]) -> None:
    """Payload refresh passes --models for single-model runs (merge mode) and omits
    it for full regens; the final command runs the payload shape tests either way.
    """
    checks = ingest.Checklist()
    ingest.run_payload_refresh(checks, model=Model.mace_mpa_0)
    *script_calls, test_call = run_cmd_calls
    assert len(script_calls) == len(ingest.PAYLOAD_SCRIPTS)
    assert any(
        "single_model_per_element_errors.py" in " ".join(cmd) for cmd in script_calls
    )
    for cmd in script_calls:
        assert cmd[2] == "--with-editable"
        assert cmd[3].startswith(".")
        assert cmd[-2:] == ("--models", "mace_mpa_0")
    kappa_call = next(
        cmd
        for cmd in script_calls
        if any("kappa_103_analysis.py" in arg for arg in cmd)
    )
    assert kappa_call[2:4] == ("--with-editable", ".[phonons]")
    assert "--extra" not in kappa_call
    assert "pytest" in test_call
    assert test_call[2:4] == ("--with-editable", ".")
    assert checks.n_failed == 0

    run_cmd_calls.clear()
    ingest.run_payload_refresh(ingest.Checklist())
    assert all("--models" not in cmd for cmd in run_cmd_calls)


def test_run_model_steps_installs_project_extras(
    run_cmd_calls: list[tuple[str, ...]],
) -> None:
    """Per-model eval subprocesses resolve project extras such as phonons."""
    checks = ingest.Checklist()
    step_cmd = "--extra phonons python scripts/evals/kappa.py"
    steps = [("Kappa metrics", True, True, step_cmd)]
    ingest.run_model_steps("evals", steps, Model.mace_mpa_0, checks, energy_only=False)

    expected = uv_cmd(
        ".[phonons]", "python", "scripts/evals/kappa.py", "--models", "mace_mpa_0"
    )
    assert run_cmd_calls == [expected]
    assert checks.n_failed == 0


def uv_cmd(project_req: str, *args: str) -> tuple[str, ...]:
    """Expected uv run command with editable project requirement."""
    return ("uv", "run", "--with-editable", project_req, *args)


@pytest.mark.parametrize(
    ("args", "expected"),
    [
        (
            '--extra phonons python -c "print(1, 2)"',
            uv_cmd(".[phonons]", "python", "-c", "print(1, 2)"),
        ),
        (
            "python child.py --extra keep --flag value",
            uv_cmd(".", "python", "child.py", "--extra", "keep", "--flag", "value"),
        ),
        (
            "--extra=phonons --extra symmetry -- python child.py --extra keep",
            uv_cmd(".[phonons,symmetry]", "python", "child.py", "--extra", "keep"),
        ),
    ],
)
def test_uv_run_args_parses_only_top_level_extras(
    args: str, expected: tuple[str, ...]
) -> None:
    """uv_run_args preserves child command args while resolving top-level extras."""
    assert ingest.uv_run_args(args) == expected


@pytest.mark.parametrize("args", ["--extra", "--extra=", '--extra ""'])
def test_uv_run_args_rejects_empty_extra(args: str) -> None:
    """Top-level --extra must have a non-empty value."""
    with pytest.raises(ValueError, match="--extra requires"):
        ingest.uv_run_args(args)


def test_map_yaml_paths() -> None:
    """--map-yaml-paths maps PR-diff YAML paths to enum names, fails on unknowns."""
    path = Model.mace_mpa_0.yaml_path.split("/models/")[-1]
    assert ingest.map_yaml_paths([f"models/{path}", f"models/{path}"]) == ["mace_mpa_0"]
    with pytest.raises(SystemExit, match="No Model enum entry"):
        ingest.map_yaml_paths(["models/new-arch/unregistered.yml"])
