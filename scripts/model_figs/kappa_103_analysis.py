"""Export per-material kappa-103 diagnostics for the /tasks/phonons page.

Writes site/src/figs/kappa-103-analysis.jsonl with, for every model with kappa_103
predictions:
- per-material SRME and scalar conductivity (vs the DFT reference values), driving a
  "kappa error vs kappa magnitude" scatter (colored by crystal system client-side)
- per-material failure flags (imaginary phonon modes, broken symmetry during
  relaxation, relaxation hit max steps), shown in the inspector scatter tooltips
- phonon frequency spectrum comparison vs DFT: per-material Wasserstein-1 distances
  plus quantile-quantile parity pairs (ML and DFT q-meshes differ - full vs
  irreducible - so frequencies are compared as weighted distributions, not pointwise)
"""

# %%
from typing import Any, Final

import numpy as np
import pandas as pd
from pymatviz.enums import Key

from matbench_discovery import figs
from matbench_discovery.cli import cli_args, is_full_model_run
from matbench_discovery.enums import DataFiles, MbdKey, Model

KAPPA_DECIMALS: Final = 4
FREQ_DECIMALS: Final = 2  # THz
# quantile levels for the plotted spectrum parity pairs (per material)
QQ_LEVELS: Final = np.linspace(0, 1, 17)
FAILURE_FLAG_COLS: Final = {
    "imag_modes": Key.has_imag_ph_modes,
    "broken_sym": "broken_symmetry",
    "max_steps": "reached_max_steps",
}


def first_scalar(value: object) -> float | None:
    """Return the first finite scalar rounded to the payload precision.

    Kappa arrays are per-temperature with 300K first; scalars pass through.
    """
    if value is None:
        return None
    arr = np.ravel(np.asarray(value, dtype=float))
    if arr.size == 0 or not np.isfinite(arr[0]):
        return None
    return round(float(arr[0]), KAPPA_DECIMALS)


def canonical_spacegroup_numbers(
    df_reference: pd.DataFrame, material_ids: list[str]
) -> list[int]:
    """Return validated canonical space groups aligned with material IDs."""
    spacegroup_numbers: list[int] = []
    for material_id in material_ids:
        value = df_reference.loc[material_id].get(Key.init_spg_num)
        try:
            numeric_value = float(value)
        except (TypeError, ValueError):
            numeric_value = float("nan")
        if (
            isinstance(value, bool | np.bool_)
            or not np.isfinite(numeric_value)
            or not numeric_value.is_integer()
            or not 1 <= numeric_value <= 230
        ):
            raise ValueError(
                f"Invalid canonical space-group number for {material_id}: {value!r}"
            )
        spacegroup_numbers.append(int(numeric_value))
    return spacegroup_numbers


def row_flag(row: pd.Series, col: str) -> bool | None:
    """A boolean failure flag from a kappa prediction row, None if absent."""
    value = row.get(col)
    if value is None or pd.isna(value):
        return None
    return bool(value)


def model_payload(
    model: Model, df_ml: pd.DataFrame, df_dft: pd.DataFrame, material_ids: list[str]
) -> dict[str, Any]:
    """Per-material diagnostics for one model, aligned to the DFT material order."""
    from matbench_discovery.metrics import phonons as phonon_metrics

    df_metrics = phonon_metrics.calc_kappa_metrics_from_dfs(df_ml.copy(), df_dft)

    # per-material columns, all aligned to material_ids (None = missing/unusable).
    # dict order defines the payload key order
    payload_cols = ("kappa_ml", "srme", "srme_censored", *FAILURE_FLAG_COLS, "freq_w1")
    cols: dict[str, list[Any]] = {col: [] for col in payload_cols}
    qq_dft: list[float] = []
    qq_ml: list[float] = []

    for mat_id in material_ids:
        if mat_id not in df_metrics.index:
            for col_list in cols.values():
                col_list.append(None)
            continue
        row_ml = df_metrics.loc[mat_id]
        row_dft = df_dft.loc[mat_id]
        cols["kappa_ml"].append(first_scalar(row_ml.get(MbdKey.kappa_tot_avg)))
        cols["srme"].append(first_scalar(row_ml.get(Key.srme)))
        cols["srme_censored"].append(row_flag(row_ml, phonon_metrics.SRME_CENSORED))
        for key, flag_col in FAILURE_FLAG_COLS.items():
            cols[key].append(row_flag(row_ml, flag_col))

        ml_freqs = phonon_metrics.phonon_spectrum_data(row_ml)
        dft_freqs = phonon_metrics.phonon_spectrum_data(row_dft)
        if ml_freqs is None or dft_freqs is None:
            cols["freq_w1"].append(None)
            continue
        ml_weights = phonon_metrics.mode_weights_for_freqs(*ml_freqs)
        dft_weights = phonon_metrics.mode_weights_for_freqs(*dft_freqs)
        # exact Wasserstein-1 distance (THz) between the ML and DFT frequency
        # spectra as weighted distributions (meshes differ: full vs irreducible, so
        # frequencies can't be paired pointwise). ML files store phono3py's BZ grid
        # whose duplicate zone-boundary points are slightly over-weighted under the
        # uniform-weight fallback -> ~0.02 THz noise floor that isn't model error
        w1_dist = phonon_metrics.phonon_spectrum_w1(row_ml, row_dft)
        cols["freq_w1"].append(first_scalar(w1_dist))
        for freqs, weights, qq_target in (
            (dft_freqs[0], dft_weights, qq_dft),
            (ml_freqs[0], ml_weights, qq_ml),
        ):
            quants = phonon_metrics.weighted_quantiles(
                freqs.ravel(), weights, QQ_LEVELS
            )
            qq_target.extend(round(float(val), FREQ_DECIMALS) for val in quants)

    finite_w1 = [val for val in cols["freq_w1"] if val is not None]
    return {
        "key": model.key,
        "label": model.label,
        **cols,
        "freq_w1_mean": round(float(np.mean(finite_w1)), KAPPA_DECIMALS)
        if finite_w1
        else None,
        "freq_pairs": {"dft": qq_dft, "ml": qq_ml},
    }


def main() -> int:
    """Build the kappa-103-analysis payload across all requested models."""
    from matbench_discovery.phonons import read_kappa_json

    df_dft = read_kappa_json(DataFiles.phonondb_pbe_103_kappa_no_nac.path)
    material_ids = [str(mat_id) for mat_id in df_dft.index]
    spg_nums = canonical_spacegroup_numbers(df_dft, material_ids)

    models: list[dict[str, Any]] = []
    for model in cli_args.models:
        try:
            kappa_path = model.kappa_103_path
        except (ValueError, FileNotFoundError):
            kappa_path = None
        if kappa_path is None:
            print(f"Skipping {model.label}: no kappa_103 predictions")
            continue
        try:
            df_ml = read_kappa_json(kappa_path)
            models.append(model_payload(model, df_ml, df_dft, material_ids))
            print(f"Processed {model.label}")
        except (ValueError, KeyError, OSError) as exc:
            print(f"Skipping {model.label}: {exc!r}")

    if not models:
        print("No models with kappa_103 predictions found")
        # on subset runs (e.g. ingesting an energy-only model) there's nothing to
        # merge, which is fine; a full run yielding no models is a config error
        return 1 if is_full_model_run() else 0

    payload = {
        "material_ids": material_ids,
        "formulas": [str(df_dft.loc[mid].get("name", "")) for mid in material_ids],
        "spg_nums": spg_nums,
        "kappa_dft": [
            first_scalar(df_dft.loc[mid].get(MbdKey.kappa_tot_avg))
            for mid in material_ids
        ],
        "models": models,
    }
    figs.write_site_payload("kappa-103-analysis", payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
