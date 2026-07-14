"""Tests for the data-only PR overlay used by CI model ingestion."""

from pathlib import Path

import pytest

import scripts.apply_pr_models_overlay as overlay


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

    assert (
        overlay.main(str(pr_root), ["models/new/missing.yml"], ["models/old/model.yml"])
        == 1
    )
    assert removed_file.is_file()

    assert (
        overlay.main(
            str(pr_root),
            [
                "models/new/model.yml",
                r"models\new\model.yml",
                "models/other/model.yml",
                "models/new/model.yml",
            ],
            ["models/old/model.yml"],
        )
        == 0
    )
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
    """A canonical-looking YAML symlink cannot cross the trust boundary."""
    trusted_root, pr_root = overlay_roots
    submitted_model = pr_root / "models/arch/model.yml"
    submitted_model.parent.mkdir(parents=True)
    try:
        submitted_model.symlink_to(pr_root / "payload.yml")
    except OSError as exc:
        pytest.skip(f"Symlinks unavailable: {exc}")
    assert overlay.main(str(pr_root), ["models/arch/model.yml"]) == 1
    assert not (trusted_root / "models").is_dir()
