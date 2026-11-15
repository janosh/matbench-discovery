"""Explicit task inapplicability must survive production roster selection."""

import gzip
import json
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
        ("not_available", True, True, True),
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
    model = Model["pet_xs_pair_bo_v1"]
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
    candidate, peer = Model["pet_xs_pair_bo_v1"], Model["mace_mp_0"]
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
