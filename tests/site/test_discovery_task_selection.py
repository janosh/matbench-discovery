"""Explicit task inapplicability must survive production roster selection."""

import gzip
import json
import os
import subprocess
import sys
from argparse import Namespace
from pathlib import Path
from types import ModuleType, SimpleNamespace
from zipfile import ZipFile

import pandas as pd
import pytest
from pymatviz.enums import Key

from matbench_discovery import cli, data, enums
from matbench_discovery.enums import MbdKey, Model
from tests.utils import import_repo_script, patch_model_metadata

if "asset_helpers" not in sys.modules:
    import_repo_script("asset_helpers", "site/scripts/asset_helpers.py")
energy_assets = import_repo_script(
    "energy_parity_task_selection", "site/scripts/generate-energy-parity-assets.py"
)


@pytest.mark.parametrize(
    ("status", "has_predictions", "active", "selected"),
    [
        ("not_applicable", True, True, False),
        ("not_applicable", False, True, False),
        ("pending", True, True, True),
        ("partial", True, True, True),
        ("complete", True, True, True),
        ("not_available", True, True, False),
        (None, True, True, True),
        (None, False, True, False),
        (None, True, False, False),
    ],
)
def test_complete_models_respects_explicit_inapplicability(
    monkeypatch: pytest.MonkeyPatch,
    status: str | None,
    has_predictions: bool,
    active: bool,
    selected: bool,
) -> None:
    """A retained artifact must not override N/A; other selection rules remain."""
    model = Model["orb_g2_groupcv_v1"]
    discovery: dict[str, object] = {"status": status} if status is not None else {}
    if has_predictions:
        discovery["pred_file"] = {"name": "retained.csv.gz"}
    patch_model_metadata(
        monkeypatch,
        model,
        model.metadata
        | {
            "lifecycle": "active" if active else "pending",
            "metrics": {"discovery": discovery},
        },
    )
    monkeypatch.setattr(cli, "cli_args", Namespace(models=[model]))

    assert cli.complete_models() == ([model] if selected else [])


@pytest.fixture
def local_parity_inputs(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> tuple[ModuleType, Model, Model, list[str], Path]:
    """Provide tiny on-disk inputs to the real loader and parity asset writer."""
    candidate, peer = Model["orb_g2_groupcv_v1"], Model["mace_mp_0"]
    for model in Model:
        patch_model_metadata(
            monkeypatch, model, model.metadata | {"lifecycle": "pending"}
        )
    monkeypatch.setattr(enums, "ROOT", str(tmp_path))
    for model in (candidate, peer):
        relative = (
            f"{model.artifact_dir}/{data.artifact_filename('2026-11-16', 'discovery')}"
        )
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({Key.mat_id: ["material-1"], "e_form_per_atom": [-1.1]}).to_csv(
            path, index=False
        )
        patch_model_metadata(
            monkeypatch,
            model,
            model.metadata
            | {
                "lifecycle": "active",
                "metrics": {"discovery": {"pred_file": {"name": relative}}},
            },
        )
    reference = pd.DataFrame(
        {
            Key.mat_id: ["material-1"],
            Key.formula: ["H"],
            Key.n_sites: [1],
            MbdKey.e_form_dft: [-1.0],
            MbdKey.each_true: [0.1],
        }
    ).set_index(str(Key.mat_id), drop=False)
    monkeypatch.setattr(data, "df_wbm", reference)
    zip_path = tmp_path / "structures.zip"
    with ZipFile(zip_path, "w") as archive:
        archive.writestr(
            "material-1.extxyz",
            '1\nLattice="5 0 0 0 5 0 0 0 5" Properties=species:S:1:pos:R:3\nH 0 0 0\n',
        )
    monkeypatch.setattr(
        energy_assets,
        "DataFiles",
        SimpleNamespace(wbm_initial_atoms=SimpleNamespace(path=str(zip_path))),
    )
    manifest = tmp_path / "manifest.json"
    argv = ["--out-dir", str(tmp_path / "parity"), "--manifest", str(manifest)]
    return energy_assets, candidate, peer, argv, manifest


@pytest.mark.parametrize(
    ("status", "has_predictions", "included"),
    [
        ("not_applicable", False, False),
        ("not_applicable", True, False),
        ("pending", True, True),
        ("partial", True, True),
        ("complete", True, True),
        (None, True, True),
    ],
)
def test_default_energy_parity_excludes_only_explicit_na(
    local_parity_inputs: tuple[ModuleType, Model, Model, list[str], Path],
    monkeypatch: pytest.MonkeyPatch,
    status: str | None,
    has_predictions: bool,
    included: bool,
) -> None:
    """Default generation skips the specialist without removing other declarations."""
    generator, candidate, peer, argv, manifest_path = local_parity_inputs
    discovery = dict(candidate.metrics["discovery"]) if has_predictions else {}
    if status is not None:
        discovery["status"] = status
    patch_model_metadata(
        monkeypatch,
        candidate,
        candidate.metadata | {"metrics": {"discovery": discovery}},
    )

    generator.main(argv)

    manifest = json.loads(manifest_path.read_text())
    expected = {peer.key, candidate.key} if included else {peer.key}
    assert set(manifest["model_assets"]) == expected
    asset = manifest["model_assets"][peer.key]["parity"]["asset"]
    with gzip.open(manifest_path.parent / "parity/assets" / asset, "rt") as file:
        assert json.load(file)["model"]["e_form_pred"] == [-1.1]


@pytest.mark.parametrize("status", [None, "pending"])
def test_energy_parity_missing_predictions_still_fail(
    local_parity_inputs: tuple[ModuleType, Model, Model, list[str], Path],
    monkeypatch: pytest.MonkeyPatch,
    status: str | None,
) -> None:
    """Missing/pending coverage does not silently excuse a broken active entry."""
    generator, candidate, _peer, argv, _manifest_path = local_parity_inputs
    discovery = {"status": status} if status is not None else {}
    patch_model_metadata(
        monkeypatch,
        candidate,
        candidate.metadata | {"metrics": {"discovery": discovery}},
    )
    with pytest.raises(ValueError, match=r"metrics\.discovery\.pred_file not found"):
        generator.main(argv)


def test_targeted_energy_na_prunes_its_previous_asset(
    local_parity_inputs: tuple[ModuleType, Model, Model, list[str], Path],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Filtering a selected N/A model still removes its stale parity record."""
    generator, candidate, peer, argv, manifest_path = local_parity_inputs
    generator.main(argv)
    before = json.loads(manifest_path.read_text())
    obsolete = before["model_assets"][candidate.key]["parity"]["asset"]
    patch_model_metadata(
        monkeypatch,
        candidate,
        candidate.metadata | {"metrics": {"discovery": {"status": "not_applicable"}}},
    )

    generator.main([*argv, "--models", candidate.name])

    after = json.loads(manifest_path.read_text())
    assert set(after["model_assets"]) == {peer.key}
    assert after["model_assets"][peer.key] == before["model_assets"][peer.key]
    assert after["base"] == before["base"]
    assert not (manifest_path.parent / "parity/assets" / obsolete).exists()


@pytest.mark.parametrize("existing_peer_manifest", [False, True])
def test_targeted_new_energy_na_needs_no_wbm_or_manifest_write(
    local_parity_inputs: tuple[ModuleType, Model, Model, list[str], Path],
    monkeypatch: pytest.MonkeyPatch,
    existing_peer_manifest: bool,
) -> None:
    """An N/A target with no old entry is a no-op even when WBM is unavailable."""
    generator, candidate, _peer, argv, manifest_path = local_parity_inputs
    patch_model_metadata(
        monkeypatch,
        candidate,
        candidate.metadata | {"metrics": {"discovery": {"status": "not_applicable"}}},
    )
    if existing_peer_manifest:
        generator.main(argv)
    before = manifest_path.read_bytes() if manifest_path.exists() else None

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("N/A-only no-op must not load data or write assets/manifests")

    for name in (
        "load_df_wbm_with_preds",
        "write_structure_shards",
        "write_json_gz",
        "write_manifest",
        "prune_unreferenced_assets",
    ):
        monkeypatch.setattr(generator, name, forbidden)
    generator.main([*argv, "--models", candidate.name])
    assert (manifest_path.read_bytes() if manifest_path.exists() else None) == before
    if not existing_peer_manifest:
        assert not (manifest_path.parent / "parity").exists()


def test_targeted_energy_na_cleanup_preserves_all_shared_data_without_wbm(
    local_parity_inputs: tuple[ModuleType, Model, Model, list[str], Path],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Remove the target offline, preserving peers, shared fields and other orphans."""
    generator, candidate, peer, argv, manifest_path = local_parity_inputs
    generator.main(argv)
    before = json.loads(manifest_path.read_text())
    before["publication_note"] = {"preserve": "future shared metadata"}
    generator.write_manifest(manifest_path, before)
    asset_dir = manifest_path.parent / "parity/assets"
    obsolete = before["model_assets"][candidate.key]["parity"]["asset"]
    orphan = asset_dir / f"{generator.ASSET_PREFIX}-unrelated-orphan.json.gz"
    orphan.write_bytes(b"unrelated prior asset")
    preserved_bytes = {
        path.name: path.read_bytes()
        for path in asset_dir.iterdir()
        if path.name != obsolete
    }
    patch_model_metadata(
        monkeypatch,
        candidate,
        candidate.metadata | {"metrics": {"discovery": {"status": "not_applicable"}}},
    )
    # Cleanup must not turn into a global lifecycle refresh either.
    patch_model_metadata(monkeypatch, peer, peer.metadata | {"lifecycle": "pending"})

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("N/A cleanup must not load WBM or regenerate base/structures")

    for name in ("load_df_wbm_with_preds", "write_structure_shards", "write_json_gz"):
        monkeypatch.setattr(generator, name, forbidden)
    generator.main([*argv, "--models", candidate.name])
    assert json.loads(manifest_path.read_text()) == before | {
        "model_assets": {peer.key: before["model_assets"][peer.key]}
    }
    assert not (asset_dir / obsolete).exists()
    assert {
        path.name: path.read_bytes() for path in asset_dir.iterdir()
    } == preserved_bytes


def test_targeted_energy_na_cleanup_rejects_changed_prefix_before_writing(
    local_parity_inputs: tuple[ModuleType, Model, Model, list[str], Path],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An offline deletion must not bypass the existing manifest-prefix guard."""
    generator, candidate, _peer, argv, manifest_path = local_parity_inputs
    generator.main(argv)
    before = manifest_path.read_bytes()
    patch_model_metadata(
        monkeypatch,
        candidate,
        candidate.metadata | {"metrics": {"discovery": {"status": "not_applicable"}}},
    )

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("Prefix validation must precede WBM or asset access")

    monkeypatch.setattr(generator, "load_df_wbm_with_preds", forbidden)
    with pytest.raises(ValueError, match="prefix changed"):
        generator.main(
            [*argv, "--models", candidate.name, "--asset-prefix", "different"]
        )
    assert manifest_path.read_bytes() == before


def test_scoped_asset_pruning_keeps_referenced_and_unrelated_files(
    tmp_path: Path,
) -> None:
    """A deletion candidate is removed only if absent from the new manifest."""
    asset_dir = tmp_path / "assets"
    asset_dir.mkdir()
    for name in (
        "parity-base.json.gz",
        "parity-target.json.gz",
        "parity-other.json.gz",
    ):
        (asset_dir / name).write_bytes(name.encode())
    manifest = {"base": {"asset": "parity-base.json.gz"}, "model_assets": {}}
    removed = energy_assets.prune_unreferenced_assets(
        asset_dir,
        "parity",
        manifest,
        target_assets={"parity-base.json.gz", "parity-target.json.gz"},
    )
    assert removed == ["parity-target.json.gz"]
    assert {path.name for path in asset_dir.iterdir()} == {
        "parity-base.json.gz",
        "parity-other.json.gz",
    }


@pytest.mark.parametrize("existing_target", [False, True])
def test_na_only_energy_command_avoids_import_time_wbm_access(
    local_parity_inputs: tuple[ModuleType, Model, Model, list[str], Path],
    existing_target: bool,
) -> None:
    """A fresh interpreter must reach both N/A paths without importing WBM data."""
    generator, candidate, peer, argv, manifest_path = local_parity_inputs
    if existing_target:
        generator.main(argv)
    script_file = generator.__file__
    assert script_file is not None
    script = Path(script_file).resolve()
    code = """
import importlib.abc
from pathlib import Path
import runpy
import sys

class BlockWbmData(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'matbench_discovery.data':
            raise RuntimeError('N/A-only command imported eager WBM data module')
        return None

script = Path(sys.argv[1])
sys.meta_path.insert(0, BlockWbmData())
sys.path.insert(0, str(script.parents[2]))
sys.path.insert(0, str(script.parent))
sys.argv = sys.argv[1:]
runpy.run_path(str(script), run_name='__main__')
"""
    cache = manifest_path.parent / "empty-cache"
    result = subprocess.run(
        [sys.executable, "-c", code, str(script), *argv, "--models", candidate.name],
        text=True,
        capture_output=True,
        check=False,
        timeout=30,
        env=os.environ | {"MBD_CACHE_DIR": str(cache), "PYTHONDONTWRITEBYTECODE": "1"},
    )
    assert result.returncode == 0, result.stderr
    assert not cache.exists()
    if existing_target:
        assert set(json.loads(manifest_path.read_text())["model_assets"]) == {peer.key}
    else:
        assert not manifest_path.exists()
