"""Generate compact assets for interactive model-page energy parity plots."""

from __future__ import annotations

import argparse
import hashlib
import math
from pathlib import Path
from typing import TYPE_CHECKING, Final
from zipfile import ZipFile

from asset_helpers import (
    clean_floats,
    clean_ints,
    compact_extxyz,
    read_manifest,
    remove_model_assets,
    resolve_models,
    retained_parity_assets,
    write_json_gz,
    write_manifest,
)
from pymatviz.enums import Key

from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.enums import DataFiles, MbdKey

if TYPE_CHECKING:
    from collections.abc import Iterable

OUT_DIR: Final = "site/static/energy-parity"
MANIFEST_PATH: Final = "site/src/lib/parity/energy-parity-manifest.json"
ASSET_PREFIX: Final = "energy-parity-wbm-v1"
LOCAL_ASSET_BASE_URL: Final = "/energy-parity/assets"
STRUCTURE_SHARD_SIZE: Final = 512
STRUCTURE_BUNDLE_SIZE: Final = 50
# round energies to shrink JSON payloads and speed up client-side parsing.
# 6 decimals = 1 ueV/atom, far finer than anything the parity plots display.
ENERGY_DECIMALS: Final = 6
# round atomic positions to 0.001 A when re-serializing structures; the header
# (Lattice, Properties, material_id, ...) and any extra columns are kept verbatim.
STRUCTURE_DECIMALS: Final = 3


def parse_args() -> argparse.Namespace:
    """Parse command-line options for generating energy parity assets."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--models",
        nargs="+",
        default=[],
        help=("Model enum names or canonical keys. Defaults to active models."),
    )
    parser.add_argument("--out-dir", default=OUT_DIR)
    parser.add_argument("--manifest", default=MANIFEST_PATH)
    parser.add_argument("--asset-prefix", default=ASSET_PREFIX)
    parser.add_argument("--local-asset-base-url", default=LOCAL_ASSET_BASE_URL)
    parser.add_argument("--limit-rows", type=int, default=None)
    parser.add_argument(
        "--structure-shard-size", type=int, default=STRUCTURE_SHARD_SIZE
    )
    parser.add_argument(
        "--structure-bundle-size", type=int, default=STRUCTURE_BUNDLE_SIZE
    )
    parser.add_argument(
        "--skip-structures",
        action="store_true",
        help="Reuse matching structure bundles; targeted runs require them.",
    )
    return parser.parse_args()


def chunked(items: list[int], chunk_size: int) -> Iterable[tuple[int, list[int]]]:
    """Yield indexed chunks from a list of integer row or shard IDs."""
    if chunk_size <= 0:
        raise ValueError(f"{chunk_size=} must be positive")
    for idx, start in enumerate(range(0, len(items), chunk_size)):
        yield idx, items[start : start + chunk_size]


def remove_stale_assets(
    asset_dir: Path, asset_prefix: str, *, keep_models: bool, keep_structures: bool
) -> None:
    """Delete old generated assets before writing a fresh manifest."""
    patterns: list[str] = []
    if not keep_models:
        patterns += [
            f"{asset_prefix}-base*.json.gz",
            f"{asset_prefix}-model-*.json.gz",
        ]
    if not keep_structures:
        patterns.append(f"{asset_prefix}-structures-*.json.gz")
    for pattern in patterns:
        for path in asset_dir.glob(pattern):
            path.unlink()


def write_structure_shards(
    material_ids: list[str],
    asset_dir: Path,
    asset_prefix: str,
    shard_size: int,
    bundle_size: int,
) -> list[dict[str, str | int]]:
    """Write compressed structure bundles grouped by shard ranges."""
    zip_path = Path(DataFiles.wbm_initial_atoms.path)
    structure_bundles: list[dict[str, str | int]] = []
    with ZipFile(zip_path) as zip_file:
        shard_idxs = list(range(math.ceil(len(material_ids) / shard_size)))
        for bundle_idx, bundle_shard_idxs in chunked(shard_idxs, bundle_size):
            shards: dict[str, dict[str, str]] = {}
            for shard_idx in bundle_shard_idxs:
                start = shard_idx * shard_size
                shards[str(shard_idx)] = {
                    material_id: compact_extxyz(
                        zip_file.read(f"{material_id}.extxyz").decode(),
                        STRUCTURE_DECIMALS,
                    )
                    for material_id in material_ids[start : start + shard_size]
                }
            asset_name = f"{asset_prefix}-structures-{bundle_idx:03d}.json.gz"
            meta = write_json_gz(asset_dir / asset_name, {"shards": shards})
            structure_bundles.append(
                {
                    **meta,
                    "start_shard": bundle_shard_idxs[0],
                    "end_shard": bundle_shard_idxs[-1],
                }
            )
    return structure_bundles


def main() -> None:
    """Generate base data, model prediction assets, structures, and manifests."""
    args = parse_args()
    out_dir = Path(args.out_dir)
    manifest_path = Path(args.manifest)
    models = resolve_models(args.models)
    df_preds = load_df_wbm_with_preds(models=models, pbar=True)
    if args.limit_rows is not None:
        df_preds = df_preds.head(args.limit_rows)

    material_ids = df_preds[Key.mat_id].astype(str).tolist()
    # hash material IDs in row order to detect stale asset manifests
    material_ids_sha256 = hashlib.sha256("\n".join(material_ids).encode()).hexdigest()
    asset_dir = out_dir / "assets"
    previous_manifest = read_manifest(manifest_path)
    target_keys = {model.key for model in models} if args.models else set()
    row_identity = (len(material_ids), material_ids_sha256)
    structure_identity = {
        "schema_version": 2,
        "row_count": row_identity[0],
        "material_ids_sha256": row_identity[1],
        "structure_shard_size": args.structure_shard_size,
        "structure_bundle_size": args.structure_bundle_size,
    }
    reused_bundles = (
        previous_manifest.get("structure_bundles") or None
        if (args.models or args.skip_structures)
        and previous_manifest is not None
        and previous_manifest.get("asset_prefix") == args.asset_prefix
        and all(
            previous_manifest.get(key) == value
            for key, value in structure_identity.items()
        )
        else None
    )
    if args.models and reused_bundles is None:
        raise ValueError("Parity structures changed; run a full refresh")
    base = {
        "material_ids": material_ids,
        "formulas": df_preds[Key.formula].astype(str).tolist(),
        "n_sites": clean_ints(df_preds[Key.n_sites]),
        "e_form_true": clean_floats(df_preds[MbdKey.e_form_dft], ENERGY_DECIMALS),
        "each_true": clean_floats(df_preds[MbdKey.each_true], ENERGY_DECIMALS),
        "structure_shard_size": args.structure_shard_size,
    }
    base_meta: dict[str, str] | None = None
    model_assets: dict[str, dict[str, str]] = {}
    if args.models:
        base_meta, model_assets = retained_parity_assets(
            previous_manifest, target_keys, base, args.asset_prefix
        )
    remove_stale_assets(
        asset_dir,
        args.asset_prefix,
        keep_models=bool(args.models),
        keep_structures=reused_bundles is not None,
    )
    remove_model_assets(asset_dir, args.asset_prefix, target_keys)
    if base_meta is None:
        base_meta = write_json_gz(asset_dir / f"{args.asset_prefix}-base.json.gz", base)

    for model in models:
        asset_name = f"{args.asset_prefix}-model-{model.key}.json.gz"
        model_payload = {
            "model_key": model.key,
            "model_label": model.label,
            "e_form_pred": clean_floats(df_preds[model.key], ENERGY_DECIMALS),
        }
        meta = write_json_gz(asset_dir / asset_name, {"model": model_payload})
        model_assets[model.key] = meta

    structure_bundles = reused_bundles or write_structure_shards(
        material_ids,
        asset_dir,
        args.asset_prefix,
        args.structure_shard_size,
        args.structure_bundle_size,
    )

    manifest = {
        **structure_identity,
        "asset_prefix": args.asset_prefix,
        "local_asset_base_url": args.local_asset_base_url.rstrip("/"),
        "base": base_meta,
        "model_assets": model_assets,
        "structure_bundles": structure_bundles,
    }
    write_manifest(manifest_path, manifest)


if __name__ == "__main__":
    main()
