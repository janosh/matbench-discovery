"""Shared helpers for generating compact gzipped site assets and manifests.

Used by the energy-parity and kappa-parity asset generators. These scripts live
in the same directory, so they import this module directly by name.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Any, TypeGuard

import numpy as np
import pandas as pd

from matbench_discovery.enums import Model
from matbench_discovery.figs import write_json_gz as write_shared_json_gz

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping

CONTENT_HASH_LENGTH = 16
# v3: per-model assets nest by kind (parity/modes) instead of a sibling mode_assets
PARITY_MANIFEST_SCHEMA_VERSION = 3


def content_addressed_name(stem: str, content_sha256: str) -> str:
    """Return an immutable asset name keyed by its uncompressed JSON content."""
    return f"{stem}-{content_sha256[:CONTENT_HASH_LENGTH]}.json.gz"


def is_content_addressed_name(stem: str, name: str) -> bool:
    """Whether an asset name ends in a lowercase hexadecimal content hash."""
    return bool(
        re.fullmatch(
            rf"{re.escape(stem)}-[0-9a-f]{{{CONTENT_HASH_LENGTH}}}\.json\.gz", name
        )
    )


def is_asset_metadata(value: object, stem: str) -> TypeGuard[dict[str, str]]:
    """Whether a manifest entry has a content-addressed name and SHA-256 digest."""
    return (
        isinstance(value, dict)
        and set(value) == {"asset", "sha256"}
        and is_content_addressed_name(stem, str(value.get("asset")))
        and isinstance(sha256 := value.get("sha256"), str)
        and re.fullmatch(r"[0-9a-f]{64}", sha256) is not None
    )


def resolve_models(model_refs: Iterable[str]) -> tuple[Model, ...]:
    """Resolve model enum names or canonical keys to members."""
    return tuple(dict.fromkeys(map(Model.from_ref, model_refs))) or Model.active()


def read_manifest(path: Path) -> dict[str, Any] | None:
    """Read an existing parity manifest when present."""
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else None


def json_content_sha256(data: Mapping[str, object]) -> str:
    """Hash the canonical uncompressed JSON used for parity asset names."""
    return hashlib.sha256(
        json.dumps(dict(data), allow_nan=False, separators=(",", ":")).encode()
    ).hexdigest()


# per-model asset kinds, mapping the manifest sub-key to the asset-name infix. Every
# model has a "parity" asset; only kappa models whose run stored the harmonic sidecar
# also have "modes"
ASSET_KINDS: dict[str, str] = {"parity": "model", "modes": "modes"}
# one manifest entry: {"asset": <content-addressed name>, "sha256": <digest>}
AssetMeta = dict[str, str]
# per-model assets keyed by model, then by kind
ModelAssets = dict[str, dict[str, AssetMeta]]


def model_asset_stem(asset_prefix: str, kind: str, model_key: str) -> str:
    """Return the content-addressed stem for one model's asset of a given kind."""
    return f"{asset_prefix}-{ASSET_KINDS[kind]}-{model_key}"


def retained_parity_assets(
    manifest: dict[str, Any] | None,
    target_keys: Iterable[str],
    base: Mapping[str, object],
    asset_prefix: str,
) -> tuple[AssetMeta, ModelAssets]:
    """Retain the immutable base and peer models from a compatible manifest.

    Each retained model maps asset kind (see ``ASSET_KINDS``) to its metadata, so a
    model's parity and modes assets cannot drift apart. Targeted and inactive models
    are dropped so the caller regenerates or prunes them.
    """
    if manifest is None:
        raise FileNotFoundError("Targeted parity refresh requires an existing manifest")
    if manifest.get("schema_version") != PARITY_MANIFEST_SCHEMA_VERSION:
        raise ValueError(
            "Parity manifest predates per-model asset kinds; run a full refresh"
        )
    if manifest.get("asset_prefix") != asset_prefix:
        raise ValueError("Parity asset prefix changed; run a full refresh")
    base_asset = manifest.get("base")
    if not is_asset_metadata(base_asset, f"{asset_prefix}-base"):
        raise ValueError("Parity base asset metadata is invalid; run a full refresh")
    expected_base_name = content_addressed_name(
        f"{asset_prefix}-base", json_content_sha256(base)
    )
    if base_asset["asset"] != expected_base_name:
        raise ValueError("Parity base content changed; run a full refresh")

    # absent (not empty) means a manifest written by a different shape; treating it as
    # empty would let a targeted refresh rewrite the manifest with only its own model
    model_assets = manifest.get("model_assets")
    if not isinstance(model_assets, dict):
        raise TypeError("Parity manifest model_assets must be an object")
    active_keys = {model.key for model in Model.active()} - set(target_keys)
    retained: ModelAssets = {}
    for model_key, kinds in model_assets.items():
        if not isinstance(model_key, str) or not isinstance(kinds, dict):
            raise TypeError("Parity model assets must map string keys to objects")
        if unknown := set(kinds) - set(ASSET_KINDS):
            raise ValueError(f"Unknown parity asset kinds {sorted(unknown)}")
        if "parity" not in kinds:
            raise ValueError(f"Parity model {model_key} has no parity asset")
        for kind, asset in kinds.items():
            if not is_asset_metadata(
                asset, model_asset_stem(asset_prefix, kind, model_key)
            ):
                raise ValueError(
                    f"Parity {kind} asset metadata is invalid; run a full refresh"
                )
        if model_key in active_keys:
            retained[model_key] = {kind: dict(asset) for kind, asset in kinds.items()}
    return dict(base_asset), retained


def prune_unreferenced_assets(
    asset_dir: Path, asset_prefix: str, manifest: Mapping[str, Any]
) -> list[str]:
    """Delete generated assets the freshly built manifest no longer references.

    Pruning after the new manifest is committed keeps every referenced asset available
    if a refresh fails. It also
    collects assets a targeted refresh cannot name, such as a model deactivated since
    the last run.
    """
    keep = {
        str(manifest["base"]["asset"]),
        *(
            str(asset["asset"])
            for kinds in manifest["model_assets"].values()
            for asset in kinds.values()
        ),
        *(str(bundle["asset"]) for bundle in manifest.get("structure_bundles", ())),
    }
    removed: list[str] = []
    for path in sorted(asset_dir.glob(f"{asset_prefix}-*.json.gz")):
        if path.name not in keep:
            path.unlink()
            removed.append(path.name)
    return removed


def clean_float(value: float | np.floating, decimals: int = 6) -> float | None:
    """Round to ``decimals``, mapping NaN/inf to None (json-safe)."""
    number = round(float(value), decimals)
    return number if math.isfinite(number) else None


def clean_floats(series: pd.Series, decimals: int = 6) -> list[float | None]:
    """Vectorized ``clean_float``; maps NaN/inf/missing to None (json-safe)."""
    nums = pd.to_numeric(series, errors="coerce").round(decimals)
    return nums.astype(object).where(np.isfinite(nums), None).tolist()


def clean_ints(series: pd.Series) -> list[int | None]:
    """Vectorized integer cleaning; maps NaN/inf/missing to None (json-safe)."""
    nums = pd.to_numeric(series, errors="coerce")
    return [int(val) if np.isfinite(val) else None for val in nums]


def write_json_gz(path: Path, data: Mapping[str, object]) -> dict[str, str]:
    """Write content-addressed gzipped JSON and return release metadata."""
    payload = dict(data)
    content_sha256 = json_content_sha256(payload)
    asset_name = content_addressed_name(
        path.name.removesuffix(".json.gz"), content_sha256
    )
    asset_path = path.with_name(asset_name)
    n_bytes = write_shared_json_gz(str(asset_path), payload)
    with open(asset_path, "rb") as file:
        sha256 = hashlib.file_digest(file, "sha256").hexdigest()
    print(f"Wrote {asset_path} ({n_bytes / 1024:.1f} KiB)")
    return {"asset": asset_name, "sha256": sha256}


def write_manifest(path: Path, manifest: Mapping[str, object]) -> None:
    """Atomically replace the manifest consumed by Python and TypeScript."""
    manifest_json = json.dumps(manifest, indent=2, sort_keys=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = tempfile.NamedTemporaryFile(  # noqa: SIM115 - closed before replace
        mode="w", encoding="utf-8", dir=path.parent, delete=False
    )
    try:
        with temporary as file:
            file.write(f"{manifest_json}\n")
        os.replace(temporary.name, path)
    finally:
        if os.path.isfile(temporary.name):
            os.remove(temporary.name)
    print(f"Wrote {path}")


def compact_extxyz(text: str, decimals: int = 3) -> str:
    """Round atomic positions to ``decimals`` places and drop column padding.

    Only the per-atom ``pos`` columns are rewritten; the header (Lattice, Properties,
    material_id, pbc, ...) and any other per-atom columns are preserved verbatim.
    """
    lines = text.splitlines()
    n_atoms = int(lines[0])
    header = lines[1]
    props = next((tok for tok in header.split() if tok.startswith("Properties=")), None)
    if props is None:
        raise ValueError(f"EXTXYZ header missing Properties= token: {header!r}")
    fields = props.split("=", 1)[1].split(":")
    col = pos_start = pos_n = 0
    for idx in range(0, len(fields), 3):
        name, n_cols = fields[idx], int(fields[idx + 2])
        if name == "pos":
            pos_start, pos_n = col, n_cols
        col += n_cols
    rows = [lines[0], header]
    for atom_line in lines[2 : 2 + n_atoms]:
        tokens = atom_line.split()
        for col_idx in range(pos_start, pos_start + pos_n):
            tokens[col_idx] = str(round(float(tokens[col_idx]), decimals))
        rows.append(" ".join(tokens))
    return "\n".join(rows) + "\n"
