"""Shared helpers for generating compact gzipped site assets and manifests.

Used by the energy-parity and kappa-parity asset generators. These scripts live
in the same directory, so they import this module directly by name.
"""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from matbench_discovery.enums import Model
from matbench_discovery.figs import write_json_gz as write_shared_json_gz

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping


def resolve_models(model_refs: Iterable[str]) -> tuple[Model, ...]:
    """Resolve model enum names, canonical keys, or key aliases to members."""
    if not model_refs:
        return Model.active()
    return tuple(map(Model.from_ref, model_refs))


def active_model_assets(
    model_assets: Mapping[str, Mapping[str, str | int]],
) -> dict[str, dict[str, str | int]]:
    """Keep only manifest assets belonging to currently active model keys."""
    active_keys = {model.key for model in Model.active()}
    return {
        model_key: dict(asset)
        for model_key, asset in model_assets.items()
        if model_key in active_keys
    }


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


def asset_safe_key(model_key: str) -> str:
    """Convert a model key into a filesystem-safe asset name segment."""
    return "".join(
        char if char.isalnum() or char in "._-" else "_" for char in model_key
    ).strip("._-")


def write_json_gz(path: Path, data: Mapping[str, object]) -> dict[str, str | int]:
    """Write deterministic gzipped JSON and return release manifest metadata."""
    n_bytes = write_shared_json_gz(str(path), dict(data))
    with open(path, "rb") as file:
        sha256 = hashlib.file_digest(file, "sha256").hexdigest()
    print(f"Wrote {path} ({n_bytes / 1024:.1f} KiB)")
    return {
        "asset": path.name,
        "bytes": n_bytes,
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
