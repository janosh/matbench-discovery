"""Tests for the data-only PR overlay used by CI model ingestion."""

import json
import re
import shutil
from pathlib import Path

import pytest

import scripts.apply_pr_models_overlay as overlay
from matbench_discovery import ROOT


@pytest.fixture
def overlay_roots(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path]:
    """Create submitted and trusted roots and enter the trusted checkout."""
    trusted_root = tmp_path / "trusted"
    trusted_root.mkdir()
    monkeypatch.chdir(trusted_root)
    return trusted_root, tmp_path / "pr"


def test_overlay_validates_then_copies_only_model_data(
    overlay_roots: tuple[Path, Path],
) -> None:
    """Atomic validation; copies model YAML only; separator duplicates collapse."""
    trusted_root, pr_root = overlay_roots
    (trusted_root / "models/old").mkdir(parents=True)
    (trusted_root / "models/untouched").mkdir()
    (pr_root / "models/new").mkdir(parents=True)
    (pr_root / "models/other").mkdir()
    (pr_root / "matbench_discovery").mkdir()
    removed_file = trusted_root / "models/old/model.yml"
    removed_file.write_text("remove: true\n")
    (trusted_root / "models/untouched/model.yml").write_text("trusted: true\n")
    (pr_root / "models/new/model.yml").write_text("submitted: true\n")
    (pr_root / "models/other/model.yml").write_text("other: true\n")
    (pr_root / "matbench_discovery/enums.py").write_text("raise RuntimeError\n")

    # every rejection bails out before the trusted checkout is touched
    removal = ["models/old/model.yml"]
    assert overlay.main(str(pr_root), []) == 1  # nothing to apply
    assert overlay.main(str(pr_root), [], ["models/gone/model.yml"]) == 1  # no such
    assert overlay.main(str(pr_root), ["models/new/missing.yml"], removal) == 1
    assert removed_file.is_file()

    yaml_paths = [
        "models/new/model.yml",
        r"models\new\model.yml",
        "models/other/model.yml",
        "models/new/model.yml",
    ]
    assert overlay.main(str(pr_root), yaml_paths, removal) == 0
    assert not removed_file.is_file()
    assert (trusted_root / "models/untouched/model.yml").is_file()
    assert (trusted_root / "models/new/model.yml").read_text() == "submitted: true\n"
    assert (trusted_root / "models/other/model.yml").read_text() == "other: true\n"
    assert not (trusted_root / "matbench_discovery").is_dir()


@pytest.mark.parametrize(
    "yaml_path",
    ["models/arch/nested/model.yml", "models/../escape.yml"],
)
def test_overlay_rejects_noncanonical_paths(
    overlay_roots: tuple[Path, Path], yaml_path: str
) -> None:
    """Only two-level model YAML paths can cross the trust boundary."""
    trusted_root, pr_root = overlay_roots
    assert overlay.main(str(pr_root), [yaml_path]) == 1
    assert not (trusted_root / "models").is_dir()


def test_overlay_rejects_symlinks(overlay_roots: tuple[Path, Path]) -> None:
    """Neither a symlinked YAML nor a symlinked models dir crosses the boundary."""
    trusted_root, pr_root = overlay_roots
    submitted_model = pr_root / "models/arch/model.yml"
    submitted_model.parent.mkdir(parents=True)
    try:
        submitted_model.symlink_to(pr_root / "payload.yml")
    except OSError as exc:
        pytest.skip(f"Symlinks unavailable: {exc}")
    assert overlay.main(str(pr_root), ["models/arch/model.yml"]) == 1
    assert not (trusted_root / "models").is_dir()

    # redirecting models/ elsewhere would otherwise smuggle in a real, non-symlink file
    (pr_root / "elsewhere/arch").mkdir(parents=True)
    (pr_root / "elsewhere/arch/model.yml").write_text("smuggled: true\n")
    shutil.rmtree(pr_root / "models")
    (pr_root / "models").symlink_to(pr_root / "elsewhere", target_is_directory=True)
    assert overlay.main(str(pr_root), ["models/arch/model.yml"]) == 1
    assert not (trusted_root / "models").is_dir()


def test_ingest_guard_exempts_only_files_the_overlay_never_copies() -> None:
    """Ingestion's blocklist exempts only files whose submitted copy is never run."""
    gh_dir = f"{ROOT}/.github"
    with open(f"{gh_dir}/workflows/update-site-figs.yml", encoding="utf-8") as file:
        workflow = file.read()
    exempt_match = re.search(r"EXEMPT_PATHS='(\[.*?\])'", workflow)
    assert exempt_match, "update-site-figs.yml no longer defines EXEMPT_PATHS"
    exempt_paths = json.loads(exempt_match[1])
    # pinned exactly: growing this set widens the ingestion trust boundary
    assert sorted(exempt_paths) == [
        "matbench_discovery/calculators.py",
        "matbench_discovery/enums.py",
    ]
    for path in exempt_paths:  # exempt only because the overlay refuses to copy them
        assert overlay.MODEL_YAML_PATTERN.fullmatch(path) is None

    # Central recipe sources stay guarded so stale committed manifests cannot ship.
    with open(f"{gh_dir}/payload-generator-paths.regex", encoding="utf-8") as file:
        generator_pattern = file.read().strip()
    for path in (
        "matbench_discovery/data.py",
        "matbench_discovery/enums.py",
        "matbench_discovery/payload_numerics.py",
        "matbench_discovery/preds/__init__.py",
    ):
        assert re.fullmatch(generator_pattern, path)
