"""Tests for the model submission-ingestion checklist without running evaluations."""

import hashlib
import json
import subprocess
from pathlib import Path

import pytest

import scripts.ingest_model as ingest
import scripts.upload_model_preds_to_figshare as upload
from matbench_discovery.enums import Model


def msgs(checks: ingest.Checklist, status: str) -> list[str]:
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


@pytest.mark.parametrize("validate_runner", [True, False])
def test_force_model_discovery_pipelines_pass_checklist(
    shared_runner_calls: list[list[str]], validate_runner: bool
) -> None:
    """Calculator-backed force models use the shared discovery runner."""
    checks = ingest.Checklist()
    assert (
        ingest.check_submission(
            Model.mace_mpa_0, checks, validate_runner=validate_runner
        )
        is False
    )
    assert not msgs(checks, ingest.FAIL)
    expected_runners = (
        {
            f"{ingest.ROOT}/models/run_discovery.py",
            f"{ingest.ROOT}/models/run_kappa.py",
            f"{ingest.ROOT}/models/run_diatomics.py",
        }
        if validate_runner
        else set()
    )
    assert {command[-4] for command in shared_runner_calls} == expected_runners
    assert all(
        command[-3:] == ["--model", Model.mace_mpa_0.name, "--dry-run"]
        for command in shared_runner_calls
    )


def test_archived_discovery_models_skip_shared_runner() -> None:
    """Archived models report why shared discovery execution is unavailable."""
    checks = ingest.Checklist()
    ingest.check_submission(Model.alignn, checks)
    assert any("discovery is archived:" in msg for msg in msgs(checks, ingest.SKIP))
    assert not any("discovery model" in msg for msg in msgs(checks, ingest.FAIL))
    assert not any(
        msg.startswith("discovery uses shared runner")
        for msg in msgs(checks, ingest.PASS)
    )


@pytest.mark.parametrize(
    ("validate_runner", "should_fail"), [(True, True), (False, False)]
)
def test_unregistered_discovery_model_validation(
    monkeypatch: pytest.MonkeyPatch, validate_runner: bool, should_fail: bool
) -> None:
    """Trusted artifact validation does not require PR calculator code."""
    monkeypatch.delitem(
        ingest.CALCULATORS,  # ty: ignore[invalid-argument-type]
        Model.mace_mpa_0.name,
    )
    checks = ingest.Checklist()
    ingest.check_submission(Model.mace_mpa_0, checks, validate_runner=validate_runner)
    failures = "\n".join(msgs(checks, ingest.FAIL))
    assert ("discovery model is not registered" in failures) is should_fail
    assert ("kappa shared runner unsupported" in failures) is should_fail


@pytest.mark.parametrize("task", ["discovery", "kappa", "diatomics"])
def test_missing_shared_runner_fails_checklist(
    monkeypatch: pytest.MonkeyPatch, task: str
) -> None:
    """Calculator-backed tasks fail validation when their shared runner is absent."""
    monkeypatch.setattr(
        ingest.os.path, "isfile", lambda path: not path.endswith(f"run_{task}.py")
    )
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
        model for model in Model if model.metadata["targets"] == "E"
    )
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
        if model_failures := msgs(checks, ingest.FAIL):
            failures[model.name] = model_failures
    assert not failures, failures


def test_malformed_file_ref_fails_checklist(monkeypatch: pytest.MonkeyPatch) -> None:
    """Malformed artifact metadata becomes a checklist failure instead of raising."""
    model = Model.mace_mpa_0
    discovery = model.metadata["metrics"]["discovery"]
    monkeypatch.setitem(discovery, "pred_file", "legacy/path.csv.gz")
    checks = ingest.Checklist()
    ingest.check_submission(model, checks, validate_runner=False)
    assert any(
        "Invalid FileRef at metrics.discovery.pred_file" in message
        for message in msgs(checks, ingest.FAIL)
    )


@pytest.mark.parametrize(
    "argv",
    [
        ["definitely-not-a-model"],
        [],
        ["--payloads-only"],
        ["mace-mpa-0", "--archive"],
        ["mace-mpa-0", "--archive-only"],
        ["mace-mpa-0", "--publish-parity"],
        ["mace-mpa-0", "--full-roster"],
    ],
)
def test_cli_rejects_invalid_args(
    monkeypatch: pytest.MonkeyPatch, argv: list[str]
) -> None:
    """Unknown, missing, or unauthenticated CLI arguments fail validation."""
    monkeypatch.delenv("FIGSHARE_TOKEN", raising=False)
    with pytest.raises(SystemExit, match="2"):
        ingest.main(argv)


def test_run_payload_refresh_scopes(run_cmd_calls: list[tuple[str, ...]]) -> None:
    """Targeted and complete refreshes pass explicit scopes and emit reports."""
    models = (Model.mace_mpa_0, Model.mace_mp_0)
    assert ingest.main([*(model.name for model in models), "--payloads-only"]) == 0
    *script_calls, report_call, test_call = run_cmd_calls
    assert len(script_calls) == len(ingest.PAYLOAD_SCRIPTS)
    assert any(
        "single_model_per_element_errors.py" in " ".join(cmd) for cmd in script_calls
    )
    assert all(
        cmd[2] == "--with-editable"
        and cmd[3].startswith(".")
        and cmd[-3:] == ("--models", *(model.name for model in models))
        for cmd in script_calls
    )
    kappa_call = next(
        cmd
        for cmd in script_calls
        if any("kappa_103_analysis.py" in arg for arg in cmd)
    )
    assert kappa_call[2:4] == ("--with-editable", ".[phonons]")
    assert "--extra" not in kappa_call
    assert "summarize_payload_changes.py" in " ".join(report_call)
    assert "pytest" in test_call
    assert test_call[2:4] == ("--with-editable", ".")
    run_cmd_calls.clear()
    ingest.run_payload_refresh(ingest.Checklist())
    *payload_calls, _report_call, _test_call = run_cmd_calls
    assert len(payload_calls) == len(ingest.PAYLOAD_SCRIPTS)
    assert all(command[-1] == "--full-roster" for command in payload_calls)


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


def test_publish_parity_assets_once(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    run_cmd_calls: list[tuple[str, ...]],
) -> None:
    """Parity publication uploads only absent immutable assets, once per type."""
    asset_bytes = b"parity asset"
    entry = {
        "asset": "asset.json.gz",
        "sha256": hashlib.sha256(asset_bytes).hexdigest(),
    }
    manifest = {"base": entry, "model_assets": {}}
    asset_paths = []
    for parity_type in ("energy", "kappa"):
        manifest_path = (
            tmp_path / f"site/src/lib/parity/{parity_type}-parity-manifest.json"
        )
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
        asset_path = (
            tmp_path / f"site/static/{parity_type}-parity/assets/{entry['asset']}"
        )
        asset_path.parent.mkdir(parents=True, exist_ok=True)
        asset_path.write_bytes(asset_bytes)
        asset_paths.append(asset_path)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(ingest, "release_asset_digests", dict)
    assert ingest.main(["--publish-parity"]) == 0
    assert len(run_cmd_calls) == 2
    assert all(
        command[:4] == ("gh", "release", "upload", "v1.0.0")
        for command in run_cmd_calls
    )
    assert all("--clobber" not in command for command in run_cmd_calls)

    for asset_path in asset_paths:
        asset_path.write_bytes(b"corrupt")
    run_cmd_calls.clear()
    assert ingest.main(["--publish-parity"]) == 1
    assert not run_cmd_calls

    published = {entry["asset"]: f"sha256:{entry['sha256']}"}
    monkeypatch.setattr(ingest, "release_asset_digests", lambda: published)
    run_cmd_calls.clear()
    assert ingest.main(["--publish-parity"]) == 0
    assert not run_cmd_calls

    published[entry["asset"]] = "sha256:conflict"
    assert ingest.main(["--publish-parity"]) == 1
    assert not run_cmd_calls


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


def test_classify_yaml_changes(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Classification targets metadata changes and rejects roster/key changes."""
    model_key = Model.mace_mpa_0.key
    active: dict[str, object] = {"model_key": model_key, "lifecycle": "active"}
    monkeypatch.setattr(ingest, "ROOT", str(tmp_path))

    def classify(
        new_files: dict[str, dict[str, object]],
        old_files: dict[str, dict[str, object]],
        removed: list[str] | None = None,
    ) -> dict[str, object]:
        """Write current files and classify them against supplied HEAD metadata."""
        for relative_path, metadata in new_files.items():
            yaml_path = tmp_path / relative_path
            yaml_path.parent.mkdir(parents=True, exist_ok=True)
            yaml_path.write_text(json.dumps(metadata))
        monkeypatch.setattr(ingest, "read_head_yaml", old_files.get)
        return ingest.classify_yaml_changes(list(new_files), removed or [])

    path = "models/edit/model.yml"
    old_geo_opt = {"metrics": {"geo_opt": {"pred_file": "old.json.gz"}}}
    assert classify({path: active}, {path: active | old_geo_opt}) == {
        "full_roster": False,
        "targeted_models": [Model.mace_mpa_0.name],
    }
    old_path, new_path = "models/old/model.yml", "models/new/model.yml"
    assert classify({new_path: active}, {old_path: active}, [old_path]) == {
        "full_roster": False,
        "targeted_models": [],
    }
    assert classify({}, {old_path: active}, [old_path])["full_roster"] is True
    assert (
        classify({path: active | {"lifecycle": "aborted"}}, {path: active})[
            "full_roster"
        ]
        is True
    )
    with pytest.raises(ValueError, match="cannot change model_key"):
        classify({path: active}, {path: active | {"model_key": "old-key"}})


@pytest.mark.parametrize(
    ("stderr", "raises"),
    [
        ("fatal: path 'new.yml' does not exist in 'HEAD'", False),
        ("fatal: path 'new.yml' exists on disk, but not in 'HEAD'", False),
        ("fatal: bad object HEAD", True),
    ],
)
def test_read_head_yaml_only_ignores_missing_paths(
    monkeypatch: pytest.MonkeyPatch, stderr: str, raises: bool
) -> None:
    """Only an absent path is treated as a newly added YAML file."""
    monkeypatch.setattr(
        ingest.subprocess,
        "run",
        lambda *_args, **_kwargs: subprocess.CompletedProcess([], 128, stderr=stderr),
    )
    if raises:
        with pytest.raises(subprocess.CalledProcessError):
            ingest.read_head_yaml("new.yml")
    else:
        assert ingest.read_head_yaml("new.yml") is None


@pytest.mark.parametrize("force_reupload", [False, True])
def test_figshare_dry_run_hashes_only_for_exact_match(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    force_reupload: bool,
) -> None:
    """Dry runs hash for exact-match reporting unless forced upload bypasses it."""

    def fake_hash(_file_path: str) -> tuple[str, int]:
        """Fail if forced dry runs perform unnecessary hashing."""
        if force_reupload:
            pytest.fail("forced dry run hashed the artifact")
        return "abc123", 11

    monkeypatch.setattr(upload.figshare, "article_exists", lambda _article_id: True)
    monkeypatch.setattr(upload.figshare, "get_existing_files", lambda _article_id: {})
    monkeypatch.setattr(upload.figshare, "get_file_hash_and_size", fake_hash)
    monkeypatch.setattr(
        upload.figshare,
        "file_exists_with_same_hash",
        lambda *_args, **_kwargs: (True, 123),
    )
    monkeypatch.setattr(
        upload, "resolve_artifact_path", lambda *_args: Model.mace_mp_0.yaml_path
    )

    upload.update_one_modeling_task_article(
        "discovery",
        [Model.mace_mp_0],
        modeling_tasks={"discovery": {"label": "Discovery"}},
        dry_run=True,
        force_reupload=force_reupload,
    )

    assert (
        "Skipped (already exists with same hash): 1" in capsys.readouterr().out
    ) is not force_reupload


@pytest.mark.parametrize(
    ("yaml_name", "artifact_rel_path", "symlink_target"),
    [
        # a submitted symlink pointing outside the model directory
        ("model.yml", "models/arch/model/leak", "secret"),
        # a real file that belongs to a sibling model in the same family
        ("model-b.yml", "models/arch/model-a/preds.csv.gz", None),
    ],
)
def test_archive_rejects_artifact_escaping_model_dir(
    yaml_name: str, artifact_rel_path: str, symlink_target: str | None, tmp_path: Path
) -> None:
    """Archive paths cannot escape the model directory via symlinks or siblings."""
    family_dir = tmp_path / "models/arch"
    artifact_path = tmp_path / artifact_rel_path
    artifact_path.parent.mkdir(parents=True)
    if symlink_target is None:
        artifact_path.touch()
    else:
        try:
            artifact_path.symlink_to(tmp_path / symlink_target)
        except OSError as exc:
            pytest.skip(f"Symlinks unavailable: {exc}")

    with pytest.raises(ValueError, match="escapes model directory"):
        upload.resolve_artifact_path(
            str(family_dir / yaml_name), artifact_rel_path, str(tmp_path)
        )
