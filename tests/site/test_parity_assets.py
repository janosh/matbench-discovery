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
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import pandas as pd
import pytest
import requests

from matbench_discovery import ROOT
from matbench_discovery.enums import Model
from tests.utils import import_repo_script

asset_helpers = import_repo_script("asset_helpers", "site/scripts/asset_helpers.py")
energy_assets = import_repo_script(
    "generate_energy_parity_assets",
    "site/scripts/generate-energy-parity-assets.py",
)
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
    assert in_manifest == expected, (
        f"missing models: {sorted(expected - in_manifest)}; "
        f"stale models: {sorted(in_manifest - expected)}. To fix, {INGEST_HINT}"
    )
    prefix = manifest["asset_prefix"]
    assert asset_helpers.is_content_addressed_name(
        f"{prefix}-base", manifest["base"]["asset"]
    )
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
    entries = [manifest["base"], *manifest["model_assets"].values()]
    entries.extend(manifest.get("structure_bundles", ()))
    # a truncated assets listing fails here (never false-passes)
    missing = {entry["asset"] for entry in entries} - published_release_assets.keys()
    assert not missing, (
        f"{kind} parity manifest references unpublished release assets: "
        f"{sorted(missing)}. To fix, {INGEST_HINT}"
    )
    for entry in entries:
        released = published_release_assets[entry["asset"]]
        expected_digest = f"sha256:{entry['sha256']}"
        assert released.get("digest") == expected_digest, (
            f"{entry['asset']}: release digest {released.get('digest')!r} != "
            f"manifest digest {expected_digest!r}"
        )

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
    assert workflow.index("generate-kappa-parity-assets.py") < workflow.index(
        "--publish-parity"
    )
    assert '"${PAYLOAD_ARGS[@]}" --payloads-only' in refresh
    with open(f"{ROOT}/.github/workflows/gh-pages.yml") as file:
        deploy = file.read()
    assert (
        ".base.asset, .model_assets[].asset, ((.structure_bundles // [])[] | .asset)"
        in deploy
    )
    assert "jq -er" in deploy
    assert 'startswith($prefix + "-")' in deploy
    assert 'mapfile -t assets < "$assets_file"' in deploy
    assert "[0-9a-f]{16}\\\\.json\\\\.gz" in deploy
    assert "-*.json.gz" not in deploy
    with open(f"{ROOT}/.github/workflows/model-pr-guard.yml") as file:
        guard = file.read()
    assert "parity/(energy|kappa)-parity-manifest" in guard


def test_targeted_parity_requires_matching_base() -> None:
    """Targeted parity keeps immutable base and peer assets from an exact manifest."""
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
    base = {"material_ids": ["material-1", "material-2"]}
    base_asset = {
        "asset": asset_helpers.content_addressed_name(
            "parity-v1-base", asset_helpers.json_content_sha256(base)
        ),
        "sha256": "d" * 64,
    }
    manifest = {
        "schema_version": 2,
        "asset_prefix": "parity-v1",
        "base": base_asset,
        "model_assets": {model_key: model_asset},
    }
    assert asset_helpers.retained_parity_assets(manifest, (), base, "parity-v1") == (
        base_asset,
        {model_key: model_asset},
    )
    with pytest.raises(FileNotFoundError, match="existing manifest"):
        asset_helpers.retained_parity_assets(None, (model_key,), base, "parity-v1")
    with pytest.raises(ValueError, match="asset prefix changed"):
        asset_helpers.retained_parity_assets(manifest, (model_key,), base, "parity-v2")
    with pytest.raises(ValueError, match="predates content-addressed assets"):
        asset_helpers.retained_parity_assets(
            manifest | {"schema_version": 1}, (model_key,), base, "parity-v1"
        )
    with pytest.raises(ValueError, match="base content changed"):
        asset_helpers.retained_parity_assets(
            manifest, (model_key,), base | {"values": [1]}, "parity-v1"
        )
    invalid_sections = {
        "base": {"base": base_asset | {"sha256": "invalid"}},
        "model": {"model_assets": {model_key: model_asset | {"sha256": "invalid"}}},
    }
    for section, update in invalid_sections.items():
        with pytest.raises(ValueError, match=f"{section} asset metadata is invalid"):
            asset_helpers.retained_parity_assets(
                manifest | update, (), base, "parity-v1"
            )


def test_structure_bundle_metadata_is_complete() -> None:
    """Structure bundle metadata requires exact count, order, names, and digests."""
    bundles = [
        {
            "asset": asset_helpers.content_addressed_name(
                f"parity-v1-structures-{bundle_idx:03d}", "a" * 64
            ),
            "sha256": "b" * 64,
        }
        for bundle_idx in range(2)
    ]
    assert energy_assets.valid_structure_bundles(bundles, "parity-v1", 2)
    invalid_bundles = [
        None,
        bundles[:1],
        bundles[::-1],
        [bundles[0] | {"sha256": "invalid"}, bundles[1]],
    ]
    assert all(
        not energy_assets.valid_structure_bundles(value, "parity-v1", 2)
        for value in invalid_bundles
    )


def test_targeted_energy_requires_matching_structure_bundles(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A targeted energy refresh cannot recreate shared structure assets."""
    df_preds = pd.DataFrame({energy_assets.Key.mat_id: ["material-1"]})
    monkeypatch.setattr(
        energy_assets, "load_df_wbm_with_preds", lambda **_kwargs: df_preds
    )
    monkeypatch.setattr(energy_assets, "read_manifest", lambda _path: {})

    with pytest.raises(ValueError, match="structures changed; run a full refresh"):
        energy_assets.main(["--models", Model.mace_mp_0.name])


@pytest.mark.parametrize(
    "parse_args", [energy_assets.parse_args, kappa_assets.parse_args]
)
def test_parity_generators_reject_empty_target_lists(
    parse_args: Callable[[list[str]], object],
) -> None:
    """An explicit empty --models cannot silently become a full refresh."""
    with pytest.raises(SystemExit):
        parse_args(["--models"])


@pytest.mark.parametrize(
    "option", ["--structure-shard-size", "--structure-bundle-size"]
)
@pytest.mark.parametrize("value", ["0", "-1"])
def test_energy_rejects_non_positive_structure_sizes(option: str, value: str) -> None:
    """Structure shard and bundle sizes must be positive."""
    with pytest.raises(SystemExit):
        energy_assets.parse_args([option, value])


@pytest.mark.parametrize("missing", ["structure", "metadata"])
def test_kappa_requires_complete_reference_data(
    monkeypatch: pytest.MonkeyPatch, missing: str
) -> None:
    """Kappa generation fails instead of emitting incomplete reference rows."""
    reference = pd.DataFrame(index=["material-1"])
    structures = {} if missing == "structure" else {"material-1": {}}
    metadata = (
        {}
        if missing == "metadata"
        else {"material-1": {"formula": "H", "n_sites": 1, "spacegroup": 1}}
    )
    monkeypatch.setattr(kappa_assets, "resolve_models", lambda _refs: ())
    monkeypatch.setattr(
        kappa_assets, "load_reference", lambda: (reference, structures, metadata)
    )
    with pytest.raises(KeyError, match="material-1"):
        kappa_assets.main([])


@pytest.mark.parametrize(
    ("resolution_error", "read_error"),
    [
        (FileNotFoundError("missing"), None),
        (ValueError("checksum"), None),
        (None, OSError("corrupt")),
    ],
)
def test_targeted_kappa_fails_on_declared_artifact_errors(
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

    monkeypatch.setattr(kappa_assets, "resolve_models", lambda _refs: (BrokenModel(),))
    reference = pd.DataFrame(
        {kappa_assets.KAPPA_TOT_AVG: [[1.0]]}, index=["material-1"]
    )
    structures = {"material-1": {"sites": []}}
    metadata = {"material-1": {"formula": "H", "n_sites": 1, "spacegroup": 1}}
    monkeypatch.setattr(
        kappa_assets, "load_reference", lambda: (reference, structures, metadata)
    )
    monkeypatch.setattr(kappa_assets, "read_manifest", lambda _path: {})
    monkeypatch.setattr(
        kappa_assets,
        "retained_parity_assets",
        lambda *_args: ({"asset": "base.json.gz", "sha256": "a" * 64}, {}),
    )
    monkeypatch.setattr(kappa_assets, "remove_model_assets", lambda *_args: None)

    def read_kappa_json(_path: str) -> pd.DataFrame:
        """Return reference data or raise the configured read failure."""
        if read_error is not None:
            raise read_error
        return reference

    monkeypatch.setattr(kappa_assets, "read_kappa_json", read_kappa_json)
    artifact_error = resolution_error or read_error
    assert artifact_error is not None
    with pytest.raises(type(artifact_error), match=str(artifact_error)):
        kappa_assets.main(["--models", BrokenModel.key])
