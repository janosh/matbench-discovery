"""Shared helpers for generating compact gzipped site assets and manifests.

Used by the energy-parity and kappa-parity asset generators. These scripts live
in the same directory, so they import this module directly by name.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

from matbench_discovery.enums import Model
from matbench_discovery.figs import write_json_gz as write_shared_json_gz

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping

CONTENT_HASH_LENGTH = 16


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


def resolve_models(model_refs: Iterable[str]) -> tuple[Model, ...]:
    """Resolve model enum names or canonical keys to members."""
    return (
        tuple(dict.fromkeys(map(Model.from_ref, model_refs)))
        if model_refs
        else Model.active()
    )


def read_manifest(path: Path) -> dict[str, Any] | None:
    """Read an existing parity manifest when present."""
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else None


def retained_model_assets(
    manifest: dict[str, Any] | None,
    target_keys: Iterable[str],
    row_identity: tuple[int, str],
    asset_prefix: str,
) -> dict[str, dict[str, str]]:
    """Retain peer assets only when they align exactly with the current base rows."""
    if manifest is None:
        raise FileNotFoundError("Targeted parity refresh requires an existing manifest")
    if manifest.get("schema_version") != 2:
        raise ValueError(
            "Parity manifest predates content-addressed assets; run a full refresh"
        )
    if (
        manifest.get("row_count"),
        manifest.get("material_ids_sha256"),
    ) != row_identity:
        raise ValueError("Parity base rows changed; run a full refresh")
    if manifest.get("asset_prefix") != asset_prefix:
        raise ValueError("Parity asset prefix changed; run a full refresh")
    model_assets = manifest.get("model_assets")
    if not isinstance(model_assets, dict):
        raise TypeError("Parity manifest model_assets must be an object")
    for model_key, asset in model_assets.items():
        if not isinstance(model_key, str) or not isinstance(asset, dict):
            raise TypeError("Parity model assets must map string keys to objects")
        if not is_content_addressed_name(
            f"{asset_prefix}-model-{model_key}", str(asset.get("asset"))
        ):
            raise ValueError(
                "Parity model asset names are not content-addressed; run a full refresh"
            )
    active_keys = {model.key for model in Model.active()} - set(target_keys)
    return {
        model_key: dict(asset)
        for model_key, asset in model_assets.items()
        if model_key in active_keys
    }


def remove_model_assets(
    asset_dir: Path, asset_prefix: str, model_keys: Iterable[str]
) -> None:
    """Remove generated files replaced by a targeted model refresh."""
    for model_key in model_keys:
        for path in asset_dir.glob(f"{asset_prefix}-model-{model_key}-*.json.gz"):
            path.unlink()


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
    content = json.dumps(payload, allow_nan=False, separators=(",", ":")).encode()
    content_sha256 = hashlib.sha256(content).hexdigest()
    suffix = ".json.gz"
    if not path.name.endswith(suffix):
        raise ValueError(f"Parity asset path must end with {suffix}: {path}")
    asset_name = content_addressed_name(path.name.removesuffix(suffix), content_sha256)
    asset_path = path.with_name(asset_name)
    n_bytes = write_shared_json_gz(str(asset_path), payload)
    with open(asset_path, "rb") as file:
        sha256 = hashlib.file_digest(file, "sha256").hexdigest()
    print(f"Wrote {asset_path} ({n_bytes / 1024:.1f} KiB)")
    return {
        "asset": asset_name,
        "sha256": sha256,
    }


def write_manifest(path: Path, manifest: Mapping[str, object]) -> None:
    """Write and format the single JSON manifest consumed by Python and TypeScript."""
    manifest_json = json.dumps(manifest, indent=2, sort_keys=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"{manifest_json}\n", encoding="utf-8")
    print(f"Wrote {path}")
    format_manifest_files([path])


def format_manifest_files(paths: Iterable[Path]) -> None:
    """Format generated manifests with the site's configured formatter."""
    site_dir = Path("site")
    vp_bin = (site_dir / "node_modules/.bin/vp").resolve()
    fmt_paths = [
        str(path.relative_to(site_dir) if path.is_relative_to(site_dir) else path)
        for path in paths
    ]
    subprocess.run(
        [str(vp_bin), "fmt", "--write", *fmt_paths],
        cwd=site_dir,
        check=True,
    )


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
