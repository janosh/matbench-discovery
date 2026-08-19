"""Guards for the energy/kappa parity plot assets on model detail pages.

Model pages fetch their energy/kappa parity plot data (too large for git) from the
v1.0.0 GitHub release. Ingesting a model (`uv run --with-editable .
scripts/ingest_model.py <model> --archive`, or the 'ingest-model' PR label which runs
it in CI) generates the assets, uploads them to the release, and commits their
metadata to the parity manifests. These tests fail fast when either step was skipped
(404s model pages in production).

model-pr-guard.yml runs this file as a required PR check (authed, so the release
test never rate-limit-skips there).
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pandas as pd
import pytest
import requests

from matbench_discovery import ROOT
from matbench_discovery.enums import Model
from tests.utils import import_repo_script

asset_helpers = import_repo_script("asset_helpers", "site/scripts/asset_helpers.py")
kappa_assets = import_repo_script(
    "generate_kappa_parity_assets", "site/scripts/generate-kappa-parity-assets.py"
)

INGEST_HINT = (
    "apply the 'ingest-model' label to the submission PR (runs the full ingest in "
    "CI) or run `uv run --with-editable . scripts/ingest_model.py <model> "
    "--archive` locally"
)


def parity_manifest(kind: str) -> dict[str, Any]:
    """Load the committed energy or kappa parity asset manifest."""
    with open(
        f"{ROOT}/site/src/lib/parity/{kind}-parity-manifest.json", encoding="utf-8"
    ) as file:
        return json.load(file)


def model_has_parity_preds(model: Model, kind: str) -> bool:
    """Whether a model has the prediction file backing its `kind` parity plot."""
    section = model.metrics.get("discovery" if kind == "energy" else "phonons") or {}
    if kind == "kappa":
        section = section.get("kappa_103") or {}
    return bool(section.get("pred_file"))


@pytest.fixture(scope="module")
def published_release_assets() -> dict[str, dict[str, Any]]:
    """Asset metadata in the GitHub release serving parity data."""
    url = "https://api.github.com/repos/janosh/matbench-discovery/releases/tags/v1.0.0"
    token = os.getenv("GH_TOKEN") or os.getenv("GITHUB_TOKEN")
    headers = {"Authorization": f"Bearer {token}"} if token else {}
    try:
        response = requests.get(url, headers=headers, timeout=30)
    except requests.RequestException as exc:  # offline local dev -> skip
        pytest.skip(f"GitHub API unreachable: {exc}")
    # rate-limited (unauthenticated runs share a per-IP quota) -> no evidence either
    # way, skip. Any other HTTP error fails: a 404 means the release itself is gone,
    # exactly the drift this test exists to catch.
    if response.status_code in (403, 429):
        pytest.skip(f"GitHub API rate-limited: {response.status_code} for {url}")
    assert response.ok, f"GitHub API error {response.status_code} for {url}"
    return {asset["name"]: asset for asset in response.json()["assets"]}


@pytest.mark.parametrize("kind", ["energy", "kappa"])
def test_parity_manifest_matches_active_models(kind: str) -> None:
    """The manifest roster exactly matches active models with predictions.

    A missing entry means a model PR was merged (or is about to be) without running
    ingest; a stale entry (e.g. a renamed or superseded model) means the generators
    were never rerun to prune it, so the site ships dead weight.
    """
    expected = {mdl.key for mdl in Model.active() if model_has_parity_preds(mdl, kind)}
    manifest = parity_manifest(kind)
    assert manifest["schema_version"] == 2
    in_manifest = set(manifest["model_assets"])
    missing = expected - in_manifest
    assert not missing, (
        f"models never ingested: {sorted(missing)}. To fix, {INGEST_HINT}"
    )
    stale = in_manifest - expected
    assert not stale, (
        f"{kind} parity manifest has entries for models that are no longer active: "
        f"{sorted(stale)}. To fix, rerun site/scripts/generate-{kind}-parity-assets.py "
        "(it prunes inactive models) and commit the refreshed manifests"
    )
    prefix = manifest["asset_prefix"]
    base = manifest["base"]
    assert asset_helpers.is_content_addressed_name(f"{prefix}-base", base["asset"])
    for model_key, asset in manifest["model_assets"].items():
        assert asset_helpers.is_content_addressed_name(
            f"{prefix}-model-{model_key}", asset["asset"]
        )
    for bundle_idx, asset in enumerate(manifest.get("structure_bundles", ())):
        assert asset_helpers.is_content_addressed_name(
            f"{prefix}-structures-{bundle_idx:03d}", asset["asset"]
        )


@pytest.mark.parametrize("kind", ["energy", "kappa"])
def test_release_has_all_parity_manifest_assets(
    kind: str, published_release_assets: dict[str, dict[str, Any]]
) -> None:
    """Release assets match manifest hashes and embed canonical model keys."""
    manifest = parity_manifest(kind)
    entries = manifest["model_assets"].values()
    bundles = manifest.get("structure_bundles", ())
    entries = (manifest["base"], *entries, *bundles)
    expected = {entry["asset"] for entry in entries}
    # a truncated assets listing fails here (never false-passes)
    missing = expected - published_release_assets.keys()
    assert not missing, (
        f"{kind} parity manifest references unpublished release assets: "
        f"{sorted(missing)}. To fix, {INGEST_HINT}"
    )
    for entry in entries:
        released = published_release_assets[entry["asset"]]
        assert released.get("digest") == f"sha256:{entry['sha256']}"

    def verify_model_asset(item: tuple[str, dict[str, Any]]) -> None:
        """Download one model asset and validate its bytes and identity."""
        model_key, entry = item
        released = published_release_assets[entry["asset"]]
        response = requests.get(released["browser_download_url"], timeout=60)
        response.raise_for_status()
        compressed = response.content
        assert hashlib.sha256(compressed).hexdigest() == entry["sha256"]
        payload = json.loads(gzip.decompress(compressed))
        assert payload["model"]["model_key"] == model_key

    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(verify_model_asset, manifest["model_assets"].items()))


def test_workflows_refresh_and_deploy_exact_parity_assets() -> None:
    """Automation refreshes, tracks, and downloads exact parity assets."""
    with open(f"{ROOT}/.github/workflows/update-site-figs.yml") as file:
        workflow = file.read()
    refresh = workflow.split("- name: Refresh payloads without secrets", 1)[1].split(
        "- name: Recheck submission head before archival", 1
    )[0]
    for kind in ("energy", "kappa"):
        assert f"generate-{kind}-parity-assets.py" in refresh
    assert 'PARITY_ARGS=(--models "${MODEL_ARGS[@]}")' in refresh
    assert refresh.count('"${PARITY_ARGS[@]}"') == 2
    assert "site/src/lib/parity/{energy,kappa}-parity-manifest.json" in workflow
    assert (
        "git add matbench_discovery/enums.py models site/src/figs site/src/lib/parity"
        in workflow
    )
    assert workflow.index("Refresh payloads without secrets") < workflow.index(
        "Archive validated model artifacts"
    )
    archive = workflow.split("- name: Archive validated model artifacts", 1)[1].split(
        "- name: Check out PR head for commit", 1
    )[0]
    assert archive.count("--publish-parity") == 1
    assert '"${MODEL_ARGS[@]}" --payloads-only' in refresh
    with open(f"{ROOT}/.github/workflows/gh-pages.yml") as file:
        deploy = file.read()
    assert ".model_assets[].asset" in deploy
    assert 'patterns+=(--pattern "$asset")' in deploy
    assert '"${patterns[@]}"' in deploy
    assert 'for asset in "${assets[@]}"' in deploy
    assert "-*.json.gz" not in deploy
    with open(f"{ROOT}/.github/workflows/model-pr-guard.yml") as file:
        guard = file.read()
    assert "parity/(energy|kappa)-parity-manifest" in guard


def test_targeted_parity_requires_matching_base_rows() -> None:
    """Targeted parity keeps peers only with an existing, row-aligned manifest."""
    assert asset_helpers.resolve_models(["mace_mp_0", "mace-mp-0"]) == (
        Model.mace_mp_0,
    )
    model_key = Model.mace_mp_0.key
    model_asset = {
        "asset": asset_helpers.content_addressed_name(
            f"parity-v1-model-{model_key}", "b" * 64
        ),
        "sha256": "a" * 64,
    }
    manifest = {
        "schema_version": 2,
        "asset_prefix": "parity-v1",
        "row_count": 2,
        "material_ids_sha256": "rows",
        "model_assets": {model_key: model_asset},
    }
    row_identity = (2, "rows")
    assert asset_helpers.retained_model_assets(
        manifest, (), row_identity, "parity-v1"
    ) == {model_key: model_asset}
    with pytest.raises(FileNotFoundError, match="existing manifest"):
        asset_helpers.retained_model_assets(
            None, (model_key,), row_identity, "parity-v1"
        )
    with pytest.raises(ValueError, match="base rows changed"):
        asset_helpers.retained_model_assets(
            manifest, (model_key,), (3, "changed"), "parity-v1"
        )
    with pytest.raises(ValueError, match="asset prefix changed"):
        asset_helpers.retained_model_assets(
            manifest, (model_key,), row_identity, "parity-v2"
        )
    with pytest.raises(ValueError, match="predates content-addressed assets"):
        asset_helpers.retained_model_assets(
            manifest | {"schema_version": 1}, (model_key,), row_identity, "parity-v1"
        )

    for script_name in ("energy", "kappa"):
        with open(
            f"{ROOT}/site/scripts/generate-{script_name}-parity-assets.py"
        ) as file:
            assert "retained_model_assets(" in file.read()


@pytest.mark.parametrize(
    ("resolution_error", "read_error"),
    [
        (FileNotFoundError("missing"), None),
        (ValueError("checksum"), None),
        (None, OSError("corrupt")),
    ],
)
def test_targeted_kappa_fails_on_declared_artifact_errors(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    resolution_error: Exception | None,
    read_error: Exception | None,
) -> None:
    """A broken declared kappa artifact cannot silently delete its manifest entry."""

    class BrokenModel:
        """Minimal model whose declared kappa artifact cannot be loaded."""

        key = "broken-model"
        label = "Broken model"

        @property
        def kappa_103_path(self) -> str:
            """Raise the configured artifact resolution failure."""
            if resolution_error is not None:
                raise resolution_error
            return "prediction.json.gz"

    monkeypatch.setattr(
        kappa_assets,
        "parse_args",
        lambda: SimpleNamespace(
            out_dir=str(tmp_path),
            manifest=str(tmp_path / "manifest.json"),
            models=[BrokenModel.key],
            asset_prefix="parity-v1",
            local_asset_base_url="/assets",
        ),
    )
    monkeypatch.setattr(kappa_assets, "resolve_models", lambda _refs: (BrokenModel(),))
    reference = pd.DataFrame(
        {kappa_assets.KAPPA_TOT_AVG: [[1.0]]}, index=["material-1"]
    )
    monkeypatch.setattr(kappa_assets, "load_reference", lambda: (reference, {}, {}))
    monkeypatch.setattr(kappa_assets, "read_manifest", lambda _path: {})
    monkeypatch.setattr(kappa_assets, "retained_model_assets", lambda *_args: {})
    monkeypatch.setattr(kappa_assets, "remove_model_assets", lambda *_args: None)

    def read_kappa_json(_path: str) -> pd.DataFrame:
        """Return reference data or raise the configured read failure."""
        if read_error is not None:
            raise read_error
        return reference

    monkeypatch.setattr(kappa_assets, "read_kappa_json", read_kappa_json)
    monkeypatch.setattr(
        kappa_assets,
        "write_json_gz",
        lambda *_args: {"asset": "base.json.gz", "sha256": "a" * 64},
    )
    monkeypatch.setattr(kappa_assets, "write_manifest", lambda *_args: None)
    artifact_error = resolution_error or read_error
    assert artifact_error is not None
    with pytest.raises(type(artifact_error), match=str(artifact_error)):
        kappa_assets.main()
