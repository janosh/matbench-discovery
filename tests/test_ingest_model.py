"""Tests for the submission-ingestion checklist (scripts/ingest_model.py).

The checklist reads model YAMLs through the Model enum metadata API (replacing the
old grep/sed bash recipe), so it can be tested against the real repo metadata
without running any evals.
"""

import subprocess

import pytest

import scripts.ingest_model as ingest
from matbench_discovery.enums import Model


def msgs(checks: "ingest.Checklist", status: str) -> list[str]:
    """Return checklist messages with the requested status."""
    return [msg for stat, msg in checks.results if stat == status]


@pytest.fixture(autouse=True)
def shared_runner_calls(monkeypatch: pytest.MonkeyPatch) -> list[list[str]]:
    """Capture shared-runner subprocesses without executing model environments."""
    calls: list[list[str]] = []

    def fake_run(command: list[str], *, check: bool) -> None:
        """Record a successful checked subprocess call."""
        assert check is True
        calls.append(command)

    monkeypatch.setattr(ingest.subprocess, "run", fake_run)
    return calls


@pytest.fixture
def run_cmd_calls(monkeypatch: pytest.MonkeyPatch) -> list[tuple[str, ...]]:
    """Patch ingest.run_cmd and return captured commands."""
    calls: list[tuple[str, ...]] = []

    def fake_run_cmd(*cmd: str) -> bool:
        calls.append(cmd)
        return True

    monkeypatch.setattr(ingest, "run_cmd", fake_run_cmd)
    return calls


def test_force_model_discovery_pipelines_pass_checklist(
    monkeypatch: pytest.MonkeyPatch, shared_runner_calls: list[list[str]]
) -> None:
    """Calculator-backed force models use the shared discovery runner."""
    original_glob = ingest.Path.glob

    def glob_with_legacy_diatomics(
        path: ingest.Path, pattern: str
    ) -> list[ingest.Path]:
        """Pretend a stale diatomics script remains beside the model YAML."""
        if pattern == "test_*_diatomics.py":
            return [path / "test_legacy_diatomics.py"]
        return list(original_glob(path, pattern))

    monkeypatch.setattr(ingest.Path, "glob", glob_with_legacy_diatomics)
    checks = ingest.Checklist()
    assert ingest.check_submission(Model.mace_mpa_0, checks) is False
    assert not msgs(checks, ingest.FAIL)
    assert not any(
        "diatomics test script found" in msg for msg in msgs(checks, ingest.PASS)
    )
    assert {command[-4] for command in shared_runner_calls} == {
        f"{ingest.ROOT}/models/run_discovery.py",
        f"{ingest.ROOT}/models/run_kappa.py",
        f"{ingest.ROOT}/models/run_diatomics.py",
    }
    assert all(
        command[-3:] == ["--model", Model.mace_mpa_0.name, "--dry-run"]
        for command in shared_runner_calls
    )


@pytest.mark.parametrize(
    "model",
    [
        Model.alignn,
        Model.equflash_29m_oam,
        Model.equiformer_v3_mp,
        Model.gnome,
    ],
)
def test_archived_discovery_models_skip_shared_runner(model: Model) -> None:
    """Archived models report why shared discovery execution is unavailable."""
    checks = ingest.Checklist()
    ingest.check_submission(model, checks)
    assert any("discovery is archived:" in msg for msg in msgs(checks, ingest.SKIP))
    assert not any("discovery model" in msg for msg in msgs(checks, ingest.FAIL))
    assert not any(
        msg.startswith("discovery uses shared runner")
        for msg in msgs(checks, ingest.PASS)
    )


def test_unregistered_discovery_model_fails_checklist(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Unclassified models cannot silently bypass shared discovery validation."""
    monkeypatch.delitem(ingest.CALCULATORS, Model.mace_mpa_0.name)
    checks = ingest.Checklist()
    ingest.check_submission(Model.mace_mpa_0, checks)
    assert any(
        "discovery model is not registered with the shared runner" in msg
        for msg in msgs(checks, ingest.FAIL)
    )


@pytest.mark.parametrize("task", ["discovery", "kappa", "diatomics"])
def test_missing_shared_runner_fails_checklist(
    monkeypatch: pytest.MonkeyPatch, task: str
) -> None:
    """Calculator-backed tasks fail validation when their shared runner is absent."""

    def fake_isfile(path: str) -> bool:
        """Pretend one shared task runner is missing."""
        return not path.endswith(f"run_{task}.py")

    monkeypatch.setattr(ingest.os.path, "isfile", fake_isfile)
    checks = ingest.Checklist()
    ingest.check_submission(Model.mace_mpa_0, checks)
    failures = msgs(checks, ingest.FAIL)
    assert any(
        f"Invalid shared {task} runner configuration: shared runner not found" in msg
        for msg in failures
    )


@pytest.mark.parametrize(
    "error",
    [
        subprocess.CalledProcessError(1, ["uv", "run"]),
        OSError("runner unavailable"),
    ],
    ids=["nonzero-exit", "os-error"],
)
def test_shared_runner_process_errors_fail_checklist(
    monkeypatch: pytest.MonkeyPatch, error: Exception
) -> None:
    """Shared-runner execution errors become checklist failures."""

    def raise_process_error(*_args: object, **_kwargs: object) -> None:
        """Raise the configured subprocess failure."""
        raise error

    monkeypatch.setattr(ingest.subprocess, "run", raise_process_error)
    checks = ingest.Checklist()
    ingest.check_submission(Model.mace_mpa_0, checks)
    failures = "\n".join(msgs(checks, ingest.FAIL))
    assert "Invalid shared discovery runner configuration:" in failures
    assert "Invalid shared diatomics runner configuration:" in failures


def test_energy_only_model_skips_force_tasks() -> None:
    """targets=E models skip geo-opt/phonons/diatomics instead of failing."""
    energy_only_model = next(
        (model for model in Model if model.metadata.get("targets") == "E"), None
    )
    if energy_only_model is None:
        pytest.skip("no energy-only model in registry")
    assert energy_only_model is not None
    checks = ingest.Checklist()
    assert ingest.check_submission(energy_only_model, checks) is True
    assert not msgs(checks, ingest.FAIL)
    skips = msgs(checks, ingest.SKIP)
    assert sum("skipped (targets=E" in msg for msg in skips) >= 4


def test_all_active_models_have_required_metadata() -> None:
    """Every active model has metadata required by evals and the site."""
    failures: dict[str, list[str]] = {}
    for model in Model.active():
        checks = ingest.Checklist()
        ingest.check_submission(model, checks)
        if hard_fails := [
            msg for msg in msgs(checks, ingest.FAIL) if "test script" not in msg
        ]:
            failures[model.name] = hard_fails
    assert not failures, failures


def test_prediction_only_kappa_model_is_explicitly_unsupported() -> None:
    """A kappa artifact without a calculator is not treated as reproducible."""
    checks = ingest.Checklist()
    ingest.check_submission(Model.matris_v050_mptrj, checks)
    assert any(
        "kappa shared runner unsupported" in msg for msg in msgs(checks, ingest.FAIL)
    )


def test_checklist_summary_counts() -> None:
    """Checklist summaries report each status count."""
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
    """Archiving requires a Figshare token."""
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
    step_cmd = "--extra phonons scripts/evals/kappa.py"
    steps = [("Kappa metrics", True, True, step_cmd)]
    ingest.run_model_steps("evals", steps, Model.mace_mpa_0, checks, energy_only=False)

    expected = uv_cmd(".[phonons]", "scripts/evals/kappa.py", "--models", "mace_mpa_0")
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
            "--extra phonons --extra symmetry python child.py --extra keep",
            uv_cmd(".[phonons,symmetry]", "python", "child.py", "--extra", "keep"),
        ),
    ],
)
def test_uv_run_args_parses_only_top_level_extras(
    args: str, expected: tuple[str, ...]
) -> None:
    """uv_run_args preserves child command args while resolving top-level extras."""
    assert ingest.uv_run_args(args) == expected


@pytest.mark.parametrize("args", ["--extra", '--extra ""'])
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
