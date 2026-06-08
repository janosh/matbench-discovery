"""Generate compact assets for interactive kappa (phonon) parity plots.

Produces a DFT-vs-ML scatter of scalar lattice thermal conductivity for the 103
phononDB-PBE materials, plus per-material structures and phonon DOS (DFT and ML)
for the click-to-inspect popups. These are compact derivatives (~57 KiB/model) of
the full kappa prediction files (12-18 MB/model), which are too large to ship to
the browser. Like the energy-parity assets, the generated .json.gz are uploaded to
the GitHub release and downloaded in CI (see .github/workflows/gh-pages.yml); only
manifest.json and the TS manifest are committed.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
from pathlib import Path
from typing import TYPE_CHECKING, Final

import ase.io
import numpy as np
import pandas as pd
from ase import Atoms
from asset_helpers import (
    asset_safe_key,
    clean_float,
    compact_extxyz,
    resolve_models,
    write_json_gz,
    write_manifest,
)
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatviz.enums import Key

from matbench_discovery.enums import DataFiles, MbdKey

if TYPE_CHECKING:
    import numpy.typing as npt

    from matbench_discovery.enums import Model

OUT_DIR: Final = "site/static/kappa-parity"
MANIFEST_TS: Final = "site/src/lib/kappa-parity-manifest.ts"
ASSET_PREFIX: Final = "kappa-parity-phonondb-v1"
LOCAL_ASSET_BASE_URL: Final = "/kappa-parity/assets"
# round conductivities to 0.0001 W/mK and structure positions to 0.001 A
KAPPA_DECIMALS: Final = 4
STRUCTURE_DECIMALS: Final = 3
# phonon DOS is precomputed as a histogram of mesh frequencies (THz) so the client
# ships a small fixed-size array instead of the full per-q-point frequency mesh.
# the client (matterviz Dos) applies Gaussian smearing, so a coarse histogram is fine
DOS_BINS: Final = 128
DOS_FREQ_DECIMALS: Final = 3
DOS_DENSITY_DECIMALS: Final = 4

KAPPA_TOT_AVG: Final = str(MbdKey.kappa_tot_avg)
PH_FREQS: Final = str(Key.ph_freqs)


def parse_args() -> argparse.Namespace:
    """Parse command-line options for generating kappa parity assets."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--models",
        nargs="*",
        default=[],
        help=(
            "Model enum names, model keys, or display labels. "
            "Defaults to active models. Models without kappa_103 predictions are "
            "skipped."
        ),
    )
    parser.add_argument("--out-dir", default=OUT_DIR)
    parser.add_argument("--manifest-ts", default=MANIFEST_TS)
    parser.add_argument("--asset-prefix", default=ASSET_PREFIX)
    parser.add_argument("--local-asset-base-url", default=LOCAL_ASSET_BASE_URL)
    return parser.parse_args()


def scalar_from_avg(value: npt.ArrayLike) -> float | None:
    """Return a precomputed directionally averaged conductivity at the first temp."""
    array = np.ravel(np.asarray(value, dtype=float))
    return clean_float(array[0], KAPPA_DECIMALS) if array.size else None


def scalar_from_rta(value: npt.ArrayLike) -> float | None:
    """Directionally average an RTA conductivity tensor at the first temperature.

    Handles Voigt 6-vectors ([xx, yy, zz, yz, xz, xy]), full 3x3 tensors, and plain
    3-vectors, each optionally prefixed by a temperature axis. When present, the
    temperature axis is always leading (component dims follow), so ``arr[0]`` selects
    the first temperature, the single point we export.
    """
    arr = np.asarray(value, dtype=float)
    if arr.size == 0:
        return None
    if arr.ndim == 3:  # (n_temp, 3, 3) -> first temperature
        arr = arr[0]
    if arr.ndim == 2 and arr.shape == (3, 3):  # full tensor -> diagonal mean
        return clean_float(float(np.trace(arr) / 3), KAPPA_DECIMALS)
    # (n_temp, 6) Voigt or (n_temp, 3) -> first temp; (6,)/(3,) -> as-is
    row = arr[0] if arr.ndim == 2 else arr
    if row.ndim != 1 or row.shape[0] < 3:
        return None
    return clean_float(float(np.mean(row[:3])), KAPPA_DECIMALS)


def row_scalar_kappa(row: pd.Series) -> float | None:
    """Scalar conductivity from a kappa row, preferring the precomputed average."""
    avg = row.get(KAPPA_TOT_AVG)
    if avg is not None and (scalar := scalar_from_avg(avg)) is not None:
        return scalar
    rta = row.get(str(MbdKey.kappa_tot_rta))
    return scalar_from_rta(rta) if rta is not None else None


def phonon_dos(
    freqs: npt.ArrayLike, weights: npt.ArrayLike | None
) -> dict[str, list[float]] | None:
    """Histogram mesh phonon frequencies (THz) into a small normalized DOS.

    Frequencies are weighted by ``weights`` only when their lengths match the
    q-point axis (irreducible mesh); otherwise a uniform full mesh is assumed.
    Densities are normalized to a peak of 1 so DFT and ML curves overlay sensibly.
    """
    freq_arr = np.asarray(freqs, dtype=float)
    if freq_arr.ndim != 2 or freq_arr.size == 0:
        return None
    flat = freq_arr.ravel()
    weight_arr = None if weights is None else np.asarray(weights, dtype=float)
    if weight_arr is not None and weight_arr.shape == (freq_arr.shape[0],):
        flat_weights = np.repeat(weight_arr, freq_arr.shape[1])
    else:
        flat_weights = None
    low, high = float(flat.min()), float(flat.max())
    if not (math.isfinite(low) and math.isfinite(high)) or high <= low:
        return None
    hist, edges = np.histogram(
        flat, bins=DOS_BINS, range=(low, high), weights=flat_weights
    )
    centers = (edges[:-1] + edges[1:]) / 2
    peak = hist.max()
    densities = hist / peak if peak > 0 else hist
    return {
        "frequencies": [round(float(freq), DOS_FREQ_DECIMALS) for freq in centers],
        "densities": [round(float(den), DOS_DENSITY_DECIMALS) for den in densities],
    }


def dos_from_row(row: pd.Series) -> dict[str, list[float]] | None:
    """Compute a phonon DOS from a kappa dataframe row if frequencies are present."""
    if PH_FREQS not in row or row[PH_FREQS] is None:
        return None
    weights = row.get("weights")
    if weights is None:
        weights = row.get(str(Key.mode_weights))
    return phonon_dos(row[PH_FREQS], weights)


def load_reference() -> tuple[
    pd.DataFrame, dict[str, str], dict[str, dict[str, str | int]]
]:
    """Load DFT kappa rows plus per-material compact structures and metadata.

    ``meta[mat_id]`` carries the chemical formula, atom count, and space group
    number (the client derives crystal system + color from the space group).
    """
    df_dft = pd.read_json(DataFiles.phonondb_pbe_103_kappa_no_nac.path)
    if "mp_id" in df_dft.columns:
        df_dft = df_dft.rename(columns={"mp_id": str(Key.mat_id)})
    df_dft = df_dft.set_index(str(Key.mat_id))

    structures: dict[str, str] = {}
    meta: dict[str, dict[str, str | int]] = {}
    for atoms in ase.io.read(DataFiles.phonondb_pbe_103_structures.path, index=":"):
        mat_id = atoms.info["material_id"]
        # keep only geometry + id; drop bulky phonon-calc metadata (fc2/q-mesh/...)
        clean = Atoms(
            symbols=atoms.get_chemical_symbols(),
            positions=atoms.positions,
            cell=atoms.cell,
            pbc=atoms.pbc,
        )
        clean.info["material_id"] = mat_id
        buffer = io.StringIO()
        ase.io.write(buffer, clean, format="extxyz")
        structures[mat_id] = compact_extxyz(buffer.getvalue(), STRUCTURE_DECIMALS)
        struct = Structure.from_ase_atoms(atoms)
        meta[mat_id] = {
            "formula": struct.formula,
            "n_sites": len(struct),
            "spacegroup": SpacegroupAnalyzer(struct).get_space_group_number(),
        }
    return df_dft, structures, meta


def main() -> None:
    """Generate base (DFT) and per-model kappa parity assets plus manifests."""
    args = parse_args()
    out_dir = Path(args.out_dir)
    manifest_ts = Path(args.manifest_ts)
    asset_dir = out_dir / "assets"
    models = resolve_models(args.models)

    df_dft, structures, meta = load_reference()
    material_ids = [str(mat_id) for mat_id in df_dft.index]
    material_ids_sha256 = hashlib.sha256("\n".join(material_ids).encode()).hexdigest()

    base = {
        "material_ids": material_ids,
        "formulas": [meta.get(mid, {}).get("formula", "") for mid in material_ids],
        "n_sites": [meta.get(mid, {}).get("n_sites") for mid in material_ids],
        "spacegroups": [meta.get(mid, {}).get("spacegroup") for mid in material_ids],
        "kappa_dft": [row_scalar_kappa(df_dft.loc[mid]) for mid in material_ids],
        "structures": {
            mid: structures[mid] for mid in material_ids if mid in structures
        },
        "dft_dos": {
            mid: dos
            for mid in material_ids
            if (dos := dos_from_row(df_dft.loc[mid])) is not None
        },
    }

    # a full run (no --models) regenerates everything; a partial run keeps existing
    # model assets so submitting a single model doesn't wipe the others
    regenerate_all = not args.models
    manifest_path = out_dir / "manifest.json"
    model_assets: dict[str, dict[str, str | int]] = {}
    if regenerate_all:
        for path in asset_dir.glob(f"{args.asset_prefix}-*.json.gz"):
            path.unlink()
    elif manifest_path.is_file():
        previous = json.loads(manifest_path.read_text(encoding="utf-8"))
        model_assets = dict(previous.get("model_assets", {}))

    base_meta = write_json_gz(asset_dir / f"{args.asset_prefix}-base.json.gz", base)

    for model in models:
        kappa_path = _kappa_path(model)
        if kappa_path is None:
            print(f"Skipping {model.label}: no kappa_103 predictions")
            continue
        try:
            df_ml = pd.read_json(kappa_path).set_index(str(Key.mat_id))
        except (ValueError, KeyError, OSError) as exc:
            print(f"Skipping {model.label}: failed reading {kappa_path}: {exc!r}")
            continue

        kappa_ml = [
            row_scalar_kappa(df_ml.loc[mid]) if mid in df_ml.index else None
            for mid in material_ids
        ]
        if not any(value is not None for value in kappa_ml):
            print(f"Skipping {model.label}: no usable kappa predictions")
            continue
        model_payload = {
            "model_key": model.key,
            "model_label": model.label,
            "kappa_ml": kappa_ml,
            "ml_dos": {
                mid: dos
                for mid in material_ids
                if mid in df_ml.index
                and (dos := dos_from_row(df_ml.loc[mid])) is not None
            },
        }
        asset_name = f"{args.asset_prefix}-model-{asset_safe_key(model.key)}.json.gz"
        model_assets[model.key] = write_json_gz(
            asset_dir / asset_name, {"model": model_payload}
        )

    manifest = {
        "schema_version": 1,
        "benchmark": "phonons",
        "dataset": "phonondb",
        "asset_prefix": args.asset_prefix,
        "local_asset_base_url": args.local_asset_base_url.rstrip("/"),
        "row_count": len(material_ids),
        "material_ids_sha256": material_ids_sha256,
        "base": base_meta,
        "model_assets": model_assets,
    }
    write_manifest(
        out_dir,
        manifest_ts,
        manifest,
        export_name="kappa_parity_manifest",
        generated_by="site/scripts/generate-kappa-parity-assets.py",
    )


def _kappa_path(model: Model) -> str | None:
    """Resolve a model's kappa_103 prediction file path, or None if unavailable."""
    try:
        path = model.kappa_103_path
    except (ValueError, FileNotFoundError):
        return None
    return path if isinstance(path, str) and Path(path).is_file() else None


if __name__ == "__main__":
    main()
