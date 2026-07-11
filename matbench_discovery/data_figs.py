"""Build model-independent data-page payloads without rendering figures."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, cast

import numpy as np
import pandas as pd
import pymatviz as pmv
from pymatgen.core import Composition
from pymatviz.enums import ElemCountMode, Key
from pymatviz.utils import spg_to_crystal_sys

from matbench_discovery import figs
from matbench_discovery.enums import MbdKey

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

    from pymatgen.analysis.phase_diagram import Entry

PLOTLY_COLORS = ("#636EFA", "#EF553B", "#00CC96")


def build_wbm_hull_dist_hist(each_true: pd.Series) -> dict[str, Any]:
    """Build the stable/unstable WBM hull-distance histogram payload."""
    mean, std = each_true.mean(), each_true.std()
    value_range = (mean - 2 * std, mean + 2 * std)
    counts, bins = np.histogram(each_true, bins=150, range=value_range)
    bins = bins[1:]  # remove left-most bin edge
    return {
        "bar_width": round(float(bins[1] - bins[0]), 6),
        "stable": {
            "x": figs.round_list(bins[bins < 0]),
            "y": counts[bins < 0].tolist(),
        },
        "unstable": {
            "x": figs.round_list(bins[bins >= 0]),
            "y": counts[bins >= 0].tolist(),
        },
        "mean": round(float(mean), 5),
        "std": round(float(std), 6),
    }


def build_wbm_e_form_hist(e_form_wbm: pd.Series) -> dict[str, Any]:
    """Build the WBM uncorrected formation-energy histogram payload."""
    counts, bins = np.histogram(e_form_wbm, bins=300, range=(-5.5, 5.5))
    return {
        "x": figs.round_list(bins[:-1]),
        "y": counts.tolist(),
        "bar_width": round(float(bins[1] - bins[0]), 6),
    }


def build_mp_elemental_ref_energies(
    elemental_ref_entries: Mapping[str, Entry],
) -> dict[str, Any]:
    """Build elemental reference energies ordered by atomic number."""
    entries = sorted(
        elemental_ref_entries.values(),
        key=lambda entry: entry.composition.elements[0].number,
    )
    return {
        "x": [float(entry.composition.elements[0].number) for entry in entries],
        "y": [round(float(entry.energy_per_atom), 2) for entry in entries],
    }


def _spacegroup_sunburst_data(
    spacegroups: Sequence[int] | pd.Series,
) -> dict[str, Any]:
    """Build Plotly-compatible flat sunburst arrays without creating a figure."""
    grouped_counts: dict[str, list[tuple[int, int]]] = {}
    for spg_num, count in pd.Series(spacegroups).value_counts().sort_index().items():
        spg_num_int = int(cast("int", spg_num))
        crystal_system = spg_to_crystal_sys(spg_num_int)
        grouped_counts.setdefault(crystal_system, []).append((spg_num_int, int(count)))

    labels: list[str] = []
    parents: list[str] = []
    values: list[int] = []
    ids: list[str] = []
    for crystal_system, system_counts in sorted(grouped_counts.items()):
        labels.append(crystal_system)
        parents.append("")
        values.append(sum(count for _spg_num, count in system_counts))
        ids.append(crystal_system)
        for spg_num, count in system_counts:
            label = str(spg_num)
            labels.append(label)
            parents.append(crystal_system)
            values.append(count)
            ids.append(f"{crystal_system}/{label}")
    return {"labels": labels, "parents": parents, "values": values, "ids": ids}


def build_spacegroup_sunbursts(
    mp_spacegroups: Sequence[int] | pd.Series,
    wbm_spacegroups: Sequence[int] | pd.Series,
) -> dict[str, Any]:
    """Build MP and WBM space-group sunburst payloads."""
    return {
        "mp": _spacegroup_sunburst_data(mp_spacegroups),
        "wbm": _spacegroup_sunburst_data(wbm_spacegroups),
    }


def build_arity_hist_payload(
    mp_formulas: pd.Series,
    mp_trj_formulas: pd.Series,
    wbm_formulas: pd.Series,
) -> dict[str, Any]:
    """Build normalized MP, MPtrj, and WBM compositional-arity histograms."""
    arity_fractions = {
        label: formulas.map(Composition).map(len).value_counts().sort_index()
        / len(formulas)
        for label, formulas in (
            ("MP", mp_formulas),
            ("MPtrj", mp_trj_formulas),
            ("WBM", wbm_formulas),
        )
    }
    df_arity = pd.DataFrame(arity_fractions).query("0 < index < 7")
    return {
        "datasets": [
            {
                "label": label,
                "color": color,
                "x": figs.round_list(df_arity.index),
                "y": figs.round_list(df_arity[label]),
            }
            for label, color in zip(df_arity.columns, PLOTLY_COLORS, strict=True)
        ]
    }


def build_mp_trj_hist_payload(df_mp_trj: pd.DataFrame) -> dict[str, Any]:
    """Build the combined MPtrj energy, force, stress, magmom, and size histograms."""
    n_sites = df_mp_trj[Key.atom_nums].map(len)
    n_sites_hist, n_sites_bins = np.histogram(
        n_sites, bins=range(1, int(n_sites.max()) + 1)
    )
    cumulative = n_sites_hist.cumsum() / n_sites_hist.sum()
    return {
        "e-form": figs.histogram(df_mp_trj[MbdKey.e_form_dft], bins=300),
        "forces": figs.histogram(
            df_mp_trj[Key.forces].explode().explode().abs(), bins=300
        ),
        "stresses": figs.histogram(df_mp_trj[Key.stress_trace], bins=300),
        "magmoms": figs.histogram(df_mp_trj[Key.magmoms].dropna().explode(), bins=300),
        "n-sites": {
            "x": figs.round_list(n_sites_bins[:-1]),
            "y": n_sites_hist.tolist(),
            "bar_width": int(n_sites_bins[1] - n_sites_bins[0]),
            "cumulative": figs.round_list(cumulative),
        },
    }


def build_route_element_counts(
    df_mp: pd.DataFrame,
    df_wbm: pd.DataFrame,
    df_mp_trj: pd.DataFrame,
) -> dict[str, pd.Series]:
    """Build every route-local periodic-table element-count series."""
    counts: dict[str, pd.Series] = {}
    for dataset, formulas in (
        ("wbm", df_wbm[Key.formula]),
        ("mp", df_mp[Key.formula]),
        ("mp-trj", df_mp_trj[Key.formula]),
    ):
        # Preserve the committed route convention: three literal "NaN" MP formulas
        # were historically parsed as missing rather than sodium nitride.
        filtered_formulas = formulas[formulas != "NaN"] if dataset == "mp" else formulas
        for count_mode in (ElemCountMode.occurrence, ElemCountMode.composition):
            series = pmv.count_elements(filtered_formulas, count_mode=count_mode)
            if dataset != "mp-trj":
                series = series.astype("Int64")
            counts[f"{dataset}-element-counts-by-{count_mode}"] = series

    steps = df_wbm.index.to_series().str.split("-").str[1].astype(int)
    if not steps.between(1, 5).all():
        raise ValueError("WBM material IDs must encode steps 1 through 5")
    for batch in range(1, 6):
        counts[f"wbm-element-counts-batch={batch}"] = pmv.count_elements(
            df_wbm.loc[steps == batch, Key.formula]
        )

    arity = df_wbm[Key.formula].map(Composition).map(len)
    for n_elements, df_group in df_wbm.groupby(arity):
        counts[f"wbm-element-counts-arity={n_elements}"] = pmv.count_elements(
            df_group[Key.formula]
        )
    return counts


def build_element_counts_payload(
    mp_occurrence_counts: pd.Series,
    wbm_unique_formulas: pd.Series,
) -> dict[str, Any]:
    """Build raw and normalized MP/WBM element-count comparison series."""
    df_counts = pd.DataFrame(index=mp_occurrence_counts.index)
    df_counts["MP"] = pd.to_numeric(mp_occurrence_counts, errors="coerce").astype(float)
    df_counts["WBM"] = pmv.count_elements(
        wbm_unique_formulas, count_mode=ElemCountMode.occurrence
    )
    min_count = 10  # only show elements with at least 10 structures
    df_counts = df_counts[df_counts.sum(axis=1) > min_count]

    payload: dict[str, Any] = {}
    for normalized in (False, True):
        scale = len(wbm_unique_formulas) / 100 if normalized else 1
        scaled_counts = df_counts / scale
        payload["normalized" if normalized else "raw"] = [
            {
                "label": dataset,
                "x": [
                    str(symbol) for symbol in scaled_counts[dataset].sort_values().index
                ],
                "y": figs.round_list(scaled_counts[dataset].sort_values()),
            }
            for dataset in ("WBM", "MP")
        ]
    return payload
