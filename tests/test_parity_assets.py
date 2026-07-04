"""Guards for the energy/kappa parity plot assets on model detail pages.

Model pages fetch their energy/kappa parity plot data (too large for git) from the
v1.0.0 GitHub release. Ingesting a model (`just ingest-model <model>`, or the
'ingest-model' PR label which runs it in CI) generates the assets, uploads them to
the release, and commits their metadata to the parity manifests. These tests fail
fast when either step was skipped (which 404s model pages in production).

model-pr-guard.yml runs this file as a required PR check (authed, so the release
test never rate-limit-skips there).
"""

from __future__ import annotations

import json
import os
from typing import Any

import pytest
import requests

from matbench_discovery import ROOT
from matbench_discovery.enums import Model

INGEST_HINT = (
    "apply the 'ingest-model' label to the submission PR (runs the full ingest in "
    "CI) or run `just ingest-model <model>` locally"
)


def parity_manifest(kind: str) -> dict[str, Any]:
    """Load the committed energy or kappa parity asset manifest."""
    with open(f"{ROOT}/site/static/{kind}-parity/manifest.json") as file:
        return json.load(file)


def model_has_parity_preds(model: Model, kind: str) -> bool:
    """Whether a model has the prediction file backing its `kind` parity plot."""
    section = (model.metrics or {}).get("discovery" if kind == "energy" else "phonons")
    if not isinstance(section, dict):  # e.g. "not available"/"not applicable"
        return False
    if kind == "kappa":
        section = section.get("kappa_103") or {}
    return bool(section.get("pred_file"))


@pytest.fixture(scope="module")
def published_release_assets() -> set[str]:
    """Asset names in the GitHub release the site build downloads parity data from."""
    url = "https://api.github.com/repos/janosh/matbench-discovery/releases/tags/v1.0.0"
    token = os.getenv("GH_TOKEN") or os.getenv("GITHUB_TOKEN")
    headers = {"Authorization": f"Bearer {token}"} if token else {}
    try:
        response = requests.get(url, headers=headers, timeout=30)
        # rate-limited (unauthenticated runs share a per-IP quota) -> no evidence
        # either way, skip. Any other HTTP error fails: a 404 here means the release
        # itself is gone, exactly the drift this test exists to catch.
        if response.status_code in (403, 429):
            pytest.skip(f"GitHub API rate-limited: {response.status_code} for {url}")
        response.raise_for_status()
    except requests.HTTPError as exc:  # e.g. release tag gone -> real failure
        pytest.fail(f"GitHub API error for {url}: {exc}")
    except requests.RequestException as exc:  # offline local dev -> skip
        pytest.skip(f"GitHub API unreachable: {exc}")
    return {asset["name"] for asset in response.json()["assets"]}


@pytest.mark.parametrize("kind", ["energy", "kappa"])
def test_parity_manifest_matches_active_models(kind: str) -> None:
    """The manifest roster exactly matches active models with predictions.

    A missing entry means a model PR was merged (or is about to be) without running
    ingest; a stale entry (e.g. a renamed or superseded model) means the generators
    were never rerun to prune it, so the site ships dead weight.
    """
    expected = {mdl.key for mdl in Model.active() if model_has_parity_preds(mdl, kind)}
    in_manifest = set(parity_manifest(kind)["model_assets"])
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


@pytest.mark.parametrize("kind", ["energy", "kappa"])
def test_release_has_all_parity_manifest_assets(
    kind: str, published_release_assets: set[str]
) -> None:
    """Every asset the parity manifest references exists in the GitHub release the
    site build downloads from. Catches ingests whose `gh release upload` step never
    completed. Skips offline; model-pr-guard.yml runs it authed as a required check.
    """
    manifest = parity_manifest(kind)
    entries = manifest["model_assets"].values()
    bundles = manifest.get("structure_bundles", ())
    expected = {entry["asset"] for entry in (manifest["base"], *entries, *bundles)}
    # gh-pages.yml only downloads assets matching <asset_prefix>-*.json.gz, so a
    # published asset with a non-matching name would still 404 in production
    prefix = manifest["asset_prefix"]
    bad_names = [
        name
        for name in sorted(expected)
        if not (name.startswith(f"{prefix}-") and name.endswith(".json.gz"))
    ]
    assert not bad_names, (
        f"{kind} assets won't match the gh-pages.yml download pattern "
        f"{prefix}-*.json.gz: {bad_names}"
    )
    # a truncated assets listing fails here (never false-passes)
    missing = expected - published_release_assets
    assert not missing, (
        f"{kind} parity manifest references unpublished release assets: "
        f"{sorted(missing)}. To fix, {INGEST_HINT}"
    )
