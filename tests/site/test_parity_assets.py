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
import subprocess
import sys
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
import requests
import yaml

from matbench_discovery import ROOT
from matbench_discovery.enums import Model
from matbench_discovery.phonons.pipeline import json_ready, write_json_records
from tests.utils import import_repo_script, make_harmonic_record

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


def model_has_phonon_sidecar(model: Model, kind: str) -> bool:
    """Whether a model's kappa run stored the harmonic phonon sidecar."""
    if kind != "kappa":
        return False
    kappa = (model.metrics.get("phonons") or {}).get("kappa_103") or {}
    return bool(kappa.get("phonon_file"))


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
    assert manifest["schema_version"] == asset_helpers.PARITY_MANIFEST_SCHEMA_VERSION
    in_manifest = set(manifest["model_assets"])
    assert in_manifest == expected, (
        f"missing models: {sorted(expected - in_manifest)}; "
        f"stale models: {sorted(in_manifest - expected)}. To fix, {INGEST_HINT}"
    )
    prefix = manifest["asset_prefix"]
    assert asset_helpers.is_content_addressed_name(
        f"{prefix}-base", manifest["base"]["asset"]
    )
    for bundle_idx, asset in enumerate(manifest.get("structure_bundles", ())):
        assert asset_helpers.is_content_addressed_name(
            f"{prefix}-structures-{bundle_idx:03d}", asset["asset"]
        )
    # modes exist exactly for models whose kappa run stored the harmonic sidecar.
    # Equality, not subset: a model declaring phonon_file without a rerun would
    # otherwise ship with the mode explorer silently missing
    with_modes = {
        model_key
        for model_key, kinds in manifest["model_assets"].items()
        if "modes" in kinds
    }
    expected_modes = {
        mdl.key for mdl in Model.active() if model_has_phonon_sidecar(mdl, kind)
    }
    assert with_modes == expected_modes, (
        f"missing mode assets: {sorted(expected_modes - with_modes)}; "
        f"stale: {sorted(with_modes - expected_modes)}. To fix, {INGEST_HINT}"
    )
    for model_key, kinds in manifest["model_assets"].items():
        assert "parity" in kinds, f"{model_key} has no parity asset"
        for asset_kind, asset in kinds.items():
            stem = asset_helpers.model_asset_stem(prefix, asset_kind, model_key)
            assert asset_helpers.is_content_addressed_name(stem, asset["asset"])


@pytest.mark.parametrize("kind", ["energy", "kappa"])
def test_release_has_all_parity_manifest_assets(
    kind: str, published_release_assets: dict[str, dict[str, Any]]
) -> None:
    """Release assets match manifest hashes and embed canonical model keys."""
    manifest = parity_manifest(kind)
    entries = [
        manifest["base"],
        *(
            asset
            for kinds in manifest["model_assets"].values()
            for asset in kinds.values()
        ),
    ]
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

    parity_entries = [
        (model_key, kinds["parity"])
        for model_key, kinds in manifest["model_assets"].items()
    ]
    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(verify_model_asset, parity_entries))


def test_workflows_refresh_and_deploy_exact_parity_assets() -> None:
    """Automation refreshes, deploys, and prunes exact parity assets."""
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
        ".base.asset, (.model_assets[] | .[].asset), "
        "((.structure_bundles // [])[] | .asset)" in deploy
    )
    assert "jq -er" in deploy
    assert 'startswith($prefix + "-")' in deploy
    assert 'mapfile -t assets < "$assets_file"' in deploy
    assert "[0-9a-f]{16}\\\\.json\\\\.gz" in deploy
    assert "-*.json.gz" not in deploy
    assert "cancel-in-progress: ${{ github.event_name == 'pull_request' }}" in deploy
    cleanup = deploy.split("prune-parity-assets:", 1)[1]
    assert "needs: deploy" in cleanup
    assert "contents: write" in cleanup
    assert "github.ref == 'refs/heads/main'" in cleanup
    assert "gh pr list --state open" in cleanup
    assert "Skipping unavailable or invalid" in cleanup
    assert "now - 3600" in cleanup
    assert 'comm -23 - "$protected"' in cleanup
    assert "gh release delete-asset v1.0.0" in cleanup
    with open(f"{ROOT}/.github/workflows/model-pr-guard.yml") as file:
        guard = file.read()
    assert "parity/(energy|kappa)-parity-manifest" in guard


@pytest.mark.skipif(
    sys.platform == "win32",
    reason="runs the ubuntu-only gh-pages cleanup step under /bin/bash with a fake gh "
    "shell script on a colon-separated PATH",
)
def test_parity_cleanup_skips_stale_pr_and_deletes_only_unprotected_asset(
    tmp_path: Path,
) -> None:
    """Cleanup tolerates an inaccessible PR and subtracts protected assets."""
    with open(f"{ROOT}/.github/workflows/gh-pages.yml", encoding="utf-8") as file:
        workflow = yaml.safe_load(file)
    cleanup = workflow["jobs"]["prune-parity-assets"]["steps"][0]["run"]
    fake_gh = tmp_path / "gh"
    fake_gh.write_text(
        r"""#!/bin/bash
set -eu
if [[ "$1 $2" == "api --method" ]]; then
  endpoint=$4
  if [[ "$endpoint" == repos/stale/repo/* ]]; then
    exit 1
  fi
  manifest=${endpoint#*contents/}
  manifest=${manifest%%\?*}
  base64 < "$REPO_ROOT/$manifest" | tr -d '\n'
elif [[ "$1 $2" == "pr list" ]]; then
  printf 'stale/repo\tdeadbeef\n'
elif [[ "$1 $2" == "release view" ]]; then
  printf '%s\n' "$PROTECTED_ASSET" "$STALE_ASSET" notes.txt
elif [[ "$1 $2" == "release delete-asset" ]]; then
  printf '%s\n' "$4" >> "$DELETED_LOG"
else
  exit 2
fi
""",
        encoding="utf-8",
    )
    fake_gh.chmod(0o755)
    stale_asset = f"energy-parity-obsolete-{'0' * 16}.json.gz"
    deleted_log = tmp_path / "deleted"
    env = os.environ | {
        "DELETED_LOG": str(deleted_log),
        "GH_REPO": "janosh/matbench-discovery",
        "GITHUB_SHA": "main-sha",
        "PATH": f"{tmp_path}:{os.environ['PATH']}",
        "PROTECTED_ASSET": parity_manifest("energy")["base"]["asset"],
        "REPO_ROOT": ROOT,
        "STALE_ASSET": stale_asset,
    }
    result = subprocess.run(
        ["/bin/bash", "-e", "-o", "pipefail", "-c", cleanup],
        cwd=ROOT,
        env=env,
        capture_output=True,
        check=False,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "Skipping unavailable or invalid stale/repo@deadbeef" in result.stderr
    assert deleted_log.read_text(encoding="utf-8").splitlines() == [stale_asset]


@pytest.fixture
def refresh_manifest() -> tuple[dict[str, Any], dict[str, Any]]:
    """Build a compatible base and model asset manifest for targeted refreshes."""
    base = {"material_ids": ["material-1", "material-2"]}
    manifest = {
        "schema_version": asset_helpers.PARITY_MANIFEST_SCHEMA_VERSION,
        "asset_prefix": "parity-v1",
        "base": {
            "asset": asset_helpers.content_addressed_name(
                "parity-v1-base", asset_helpers.json_content_sha256(base)
            ),
            "sha256": "d" * 64,
        },
        "model_assets": {
            model_key: {
                kind: {
                    "asset": asset_helpers.content_addressed_name(
                        asset_helpers.model_asset_stem("parity-v1", kind, model_key),
                        "b" * 64,
                    ),
                    "sha256": "a" * 64,
                }
                for kind in ("parity", "modes")
            }
            for model_key in (Model.mace_mp_0.key, Model.chgnet_0_3_0.key, "retired")
        },
    }
    return manifest, base


def test_targeted_parity_requires_matching_base(
    refresh_manifest: tuple[dict[str, Any], dict[str, Any]],
) -> None:
    """Targeted refresh rejects missing, stale, or malformed manifest entries."""
    manifest, base = refresh_manifest
    model_key = Model.mace_mp_0.key
    kinds = manifest["model_assets"][model_key]
    with pytest.raises(FileNotFoundError, match="existing manifest"):
        asset_helpers.retained_parity_assets(None, (), base, "parity-v1")
    with pytest.raises(ValueError, match="base content changed"):
        asset_helpers.retained_parity_assets(
            manifest, (), base | {"values": [1]}, "parity-v1"
        )
    for update, message in (
        ({"asset_prefix": "parity-v2"}, "asset prefix changed"),
        ({"schema_version": 2}, "predates per-model asset kinds"),
        ({"base": manifest["base"] | {"sha256": "invalid"}}, "base asset metadata"),
        ({"model_assets": {model_key: {"modes": kinds["modes"]}}}, "has no parity"),
        ({"model_assets": {model_key: kinds | {"bogus": kinds["parity"]}}}, "Unknown"),
        (
            {
                "model_assets": {
                    model_key: {"parity": kinds["parity"] | {"sha256": "invalid"}}
                }
            },
            "parity asset metadata",
        ),
    ):
        with pytest.raises(ValueError, match=message):
            asset_helpers.retained_parity_assets(
                manifest | update, (), base, "parity-v1"
            )
    # An absent model_assets must fail: returning {} would erase peer models.
    manifest.pop("model_assets")
    with pytest.raises(TypeError, match="model_assets must be an object"):
        asset_helpers.retained_parity_assets(manifest, (), base, "parity-v1")


@pytest.mark.parametrize("targeted", [False, True])
def test_retained_parity_assets_drops_targeted_and_inactive_models(
    refresh_manifest: tuple[dict[str, Any], dict[str, Any]],
    targeted: bool,
) -> None:
    """Retain immutable base and active peers, including their complete asset kinds."""
    assert asset_helpers.resolve_models(["mace_mp_0", "mace-mp-0"]) == (
        Model.mace_mp_0,
    )
    manifest, base = refresh_manifest
    target_keys = (Model.mace_mp_0.key,) if targeted else ()
    expected = {
        key: kinds
        for key, kinds in manifest["model_assets"].items()
        if key != "retired" and key not in target_keys
    }
    assert asset_helpers.retained_parity_assets(
        manifest, target_keys, base, "parity-v1"
    ) == (
        manifest["base"],
        expected,
    )


@pytest.mark.parametrize("replace_fails", [False, True])
def test_prune_unreferenced_assets_keeps_exactly_the_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, replace_fails: bool
) -> None:
    """Pruning after the write deletes only files the new manifest does not name."""
    model_key = Model.mace_mp_0.key
    # a peer whose key extends the target key must survive: name-prefix matching alone
    # would delete `<key>-variant` when refreshing `<key>`
    extended_key = f"{model_key}-variant"
    kept = {
        "parity-v1-base-" + "0" * 16 + ".json.gz",
        f"parity-v1-model-{model_key}-{'0' * 16}.json.gz",
        f"parity-v1-modes-{model_key}-{'0' * 16}.json.gz",
        f"parity-v1-model-{extended_key}-{'0' * 16}.json.gz",
        "parity-v1-structures-000-" + "0" * 16 + ".json.gz",
    }
    stale = {
        f"parity-v1-model-{model_key}-{'1' * 16}.json.gz",  # superseded content hash
        "parity-v1-model-retired-model-" + "0" * 16 + ".json.gz",  # deactivated model
        f"parity-v1-modes-{extended_key}-{'0' * 16}.json.gz",  # dropped sidecar
    }
    for name in kept | stale:
        (tmp_path / name).touch()
    (tmp_path / "unrelated-prefix-asset.json.gz").touch()

    manifest = {
        "base": {"asset": "parity-v1-base-" + "0" * 16 + ".json.gz"},
        "model_assets": {
            model_key: {
                "parity": {"asset": f"parity-v1-model-{model_key}-{'0' * 16}.json.gz"},
                "modes": {"asset": f"parity-v1-modes-{model_key}-{'0' * 16}.json.gz"},
            },
            extended_key: {
                "parity": {
                    "asset": f"parity-v1-model-{extended_key}-{'0' * 16}.json.gz"
                }
            },
        },
        "structure_bundles": [
            {"asset": "parity-v1-structures-000-" + "0" * 16 + ".json.gz"}
        ],
    }
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text('{"previous": true}\n')
    if replace_fails:

        def fail_replace(*_args: object) -> None:
            """Simulate a failed manifest replacement."""
            raise OSError("replace failed")

        monkeypatch.setattr(asset_helpers.os, "replace", fail_replace)
        with pytest.raises(OSError, match="replace failed"):
            asset_helpers.write_manifest(manifest_path, manifest)
        assert manifest_path.read_text() == '{"previous": true}\n'
        assert {path.name for path in tmp_path.iterdir()} == kept | stale | {
            "manifest.json",
            "unrelated-prefix-asset.json.gz",
        }
        return
    asset_helpers.write_manifest(manifest_path, manifest)
    assert json.loads(manifest_path.read_text()) == manifest
    removed = asset_helpers.prune_unreferenced_assets(tmp_path, "parity-v1", manifest)
    assert set(removed) == stale
    # files under another prefix are never touched
    assert {path.name for path in tmp_path.iterdir()} == kept | {
        "manifest.json",
        "unrelated-prefix-asset.json.gz",
    }


def test_compact_phonon_modes_keeps_eigenvectors_at_labeled_points_only() -> None:
    """The site modes asset drops eigenvectors between high-symmetry points."""
    record = make_harmonic_record()
    compact = kappa_assets.compact_phonon_modes(record)
    assert set(compact) == {
        "primitive",
        "segments",
        "q_points",
        "distances",
        "frequencies",
        "eigenvectors",
    }
    assert compact["primitive"]["symbols"] == ["Na", "Cl"]
    assert compact["primitive"]["lattice"][0][0] == round(3.123456789, 4)
    assert compact["segments"] == record["band_path"]["segments"]
    # every path point keeps frequencies but eigenvectors only survive at endpoints
    assert len(compact["frequencies"]) == 4
    assert set(compact["eigenvectors"]) == {"0", "1", "2", "3"}
    np.testing.assert_allclose(
        compact["frequencies"], record["band_path"]["frequencies"], atol=5e-5, rtol=0
    )
    for q_idx, vectors in compact["eigenvectors"].items():
        np.testing.assert_allclose(
            vectors, record["band_path"]["eigenvectors"][int(q_idx)], atol=5e-5, rtol=0
        )
    # extend the path so an interior point exists and check it is dropped
    record["band_path"]["segments"][1]["end_index"] = 3
    record["band_path"]["segments"][1]["start_index"] = 1
    record["band_path"]["segments"][0]["end_index"] = 1
    compact = kappa_assets.compact_phonon_modes(record)
    assert set(compact["eigenvectors"]) == {"0", "1", "3"}
    # JSON must round-trip (no numpy scalars/arrays left behind)
    json.dumps(compact)

    # a wrong per-q-point layout (3 atoms instead of 2) and a q-count mismatch both fail
    record["band_path"]["eigenvectors"] = np.zeros((4, 6, 3, 3, 2))
    with pytest.raises(ValueError, match="have shape"):
        kappa_assets.compact_phonon_modes(record)
    record["band_path"]["eigenvectors"] = np.zeros((3, 6, 2, 3, 2))
    with pytest.raises(ValueError, match="band-path arrays disagree"):
        kappa_assets.compact_phonon_modes(record)
    record["band_path"]["segments"] = []
    with pytest.raises(ValueError, match="no band-path segments"):
        kappa_assets.compact_phonon_modes(record)


# the client zips these arrays by index without rechecking, so a length mismatch here
# truncates the band path or yields undefined displacements instead of failing loudly
@pytest.mark.parametrize(
    ("mutate", "expected"),
    [
        (
            lambda rec: rec["band_path"].__setitem__("q_points", np.zeros((2, 3))),
            "band-path arrays disagree",
        ),
        (
            lambda rec: rec["band_path"].__setitem__("distances", np.zeros(7)),
            "band-path arrays disagree",
        ),
        (
            lambda rec: rec["primitive"].__setitem__("masses", [1.0]),
            "primitive arrays disagree",
        ),
        (
            lambda rec: rec["primitive"].__setitem__("frac_coords", np.zeros((1, 3))),
            "primitive arrays disagree",
        ),
        (
            lambda rec: rec["band_path"].__setitem__("frequencies", np.zeros((4, 5))),
            "expected 6 bands for 2 atoms",
        ),
        (
            lambda rec: rec["band_path"]["segments"][1].__setitem__("end_index", 99),
            "outside 0..3",
        ),
        (
            lambda rec: rec["band_path"]["segments"][0].__setitem__("start_index", -1),
            "outside 0..3",
        ),
    ],
)
def test_compact_phonon_modes_rejects_inconsistent_shapes(
    mutate: Callable[[dict[str, Any]], None], expected: str
) -> None:
    """Every array the client indexes is length-checked against the frequency grid."""
    record = make_harmonic_record()
    mutate(record)
    with pytest.raises(ValueError, match=expected):
        kappa_assets.compact_phonon_modes(record)


def test_load_phonon_modes(tmp_path: Path) -> None:
    """A merged phonon sidecar loads into compact per-material modes keyed by ID."""
    path = tmp_path / "model-phonons.json.gz"
    with gzip.open(path, mode="wt", encoding="utf-8") as file:
        write_json_records(file, [make_harmonic_record()])
    modes = kappa_assets.load_phonon_modes(str(path))
    assert list(modes) == ["material-1"]
    assert modes["material-1"]["segments"][0]["start_label"] == "GAMMA"
    failed = make_harmonic_record(material_id="material-2")
    del failed["band_path"]
    failed["errors"] = {"band_path": {"message": "seekpath failed"}}
    with gzip.open(path, mode="wt", encoding="utf-8") as file:
        write_json_records(file, [make_harmonic_record(), failed])
    assert list(kappa_assets.load_phonon_modes(str(path))) == ["material-1"]
    with gzip.open(path, mode="wt", encoding="utf-8") as file:
        write_json_records(file, [failed, failed])
    with pytest.raises(ValueError, match="duplicate material material-2"):
        kappa_assets.load_phonon_modes(str(path))
    with gzip.open(path, mode="wt", encoding="utf-8") as file:
        write_json_records(file, [failed])
    assert kappa_assets.load_phonon_modes(str(path)) == {}
    with gzip.open(path, mode="wt", encoding="utf-8") as file:
        write_json_records(file, [])
    with pytest.raises(ValueError, match="has no records"):
        kappa_assets.load_phonon_modes(str(path))
    del failed["errors"]
    with gzip.open(path, mode="wt", encoding="utf-8") as file:
        write_json_records(file, [failed])
    with pytest.raises(KeyError, match="band_path"):
        kappa_assets.load_phonon_modes(str(path))


def test_load_phonon_modes_is_order_and_format_independent(tmp_path: Path) -> None:
    """Record order and JSON layout must not reach the content-addressed asset name."""
    records = [make_harmonic_record(material_id=f"material-{idx}") for idx in (1, 2, 3)]

    def write(name: str, payload: list[dict[str, Any]], *, plain: bool) -> str:
        path = tmp_path / name
        with gzip.open(path, mode="wt", encoding="utf-8") as file:
            if plain:  # what a contributor gets from a bare json.dump
                json.dump([json_ready(rec) for rec in payload], file, indent=2)
            else:
                write_json_records(file, payload)
        return str(path)

    framed = kappa_assets.load_phonon_modes(write("a.json.gz", records, plain=False))
    reversed_order = kappa_assets.load_phonon_modes(
        write("b.json.gz", records[::-1], plain=False)
    )
    indented = kappa_assets.load_phonon_modes(write("c.json.gz", records, plain=True))
    assert list(framed) == ["material-1", "material-2", "material-3"]
    # identical content must serialize identically regardless of input order/layout
    assert framed == reversed_order == indented

    duplicated = write("d.json.gz", [*records, records[0]], plain=False)
    with pytest.raises(ValueError, match="duplicate material material-1"):
        kappa_assets.load_phonon_modes(duplicated)
    not_an_array = tmp_path / "e.json.gz"
    with gzip.open(not_an_array, mode="wt", encoding="utf-8") as file:
        json.dump({"material-1": "a mapping, not the expected array"}, file)
    with pytest.raises(TypeError, match="must hold a JSON array"):
        kappa_assets.load_phonon_modes(str(not_an_array))


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
    monkeypatch.setattr(
        kappa_assets, "read_manifest", lambda _path: {"model_assets": {}}
    )
    monkeypatch.setattr(
        kappa_assets,
        "retained_parity_assets",
        lambda *_args: ({"asset": "base.json.gz", "sha256": "a" * 64}, {}),
    )
    monkeypatch.setattr(
        kappa_assets, "prune_unreferenced_assets", lambda *_args, **_kwargs: []
    )

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
