"""Tests for model-independent data-page payload builders."""

import gzip
import json

import pandas as pd
import pytest
from pymatviz.enums import Key

from matbench_discovery import ROOT
from matbench_discovery.data_figs import (
    SERIES_COLORS,
    build_arity_hist_payload,
    build_element_counts_payload,
    build_mp_elemental_ref_energies,
    build_mp_trj_hist_payload,
    build_route_element_counts,
    build_spacegroup_sunbursts,
    build_wbm_e_form_hist,
    build_wbm_hull_dist_hist,
)
from matbench_discovery.energy import mp_elem_ref_entries
from matbench_discovery.enums import DataFiles, MbdKey


def test_build_wbm_hull_dist_hist() -> None:
    """Hull-distance payload separates stable and unstable bins."""
    payload = build_wbm_hull_dist_hist(pd.Series([-1.0, -0.2, 0.1, 0.8]))
    assert payload["bar_width"] > 0
    assert payload["mean"] == pytest.approx(-0.075)
    assert len(payload["stable"]["x"]) == len(payload["stable"]["y"])
    assert len(payload["unstable"]["x"]) == len(payload["unstable"]["y"])


def test_build_wbm_e_form_hist() -> None:
    """Formation-energy payload uses fixed bins and excludes out-of-range values."""
    payload = build_wbm_e_form_hist(pd.Series([-6.0, -1.0, 0.0, 1.0, 6.0]))
    assert len(payload["x"]) == len(payload["y"]) == 300
    assert sum(payload["y"]) == 3
    assert payload["bar_width"] == 0.036667


def test_build_spacegroup_sunbursts() -> None:
    """Sunburst builders emit flat child and crystal-system arrays."""
    payload = build_spacegroup_sunbursts(
        pd.Series([1, 1, 2, 225]), pd.Series([225, 225, 1])
    )
    assert payload["mp"] == {
        "labels": ["cubic", "225", "triclinic", "1", "2"],
        "parents": ["", "cubic", "", "triclinic", "triclinic"],
        "values": [1, 1, 3, 2, 1],
        "ids": [
            "cubic",
            "cubic/225",
            "triclinic",
            "triclinic/1",
            "triclinic/2",
        ],
    }


def test_build_arity_hist_payload() -> None:
    """Arity payload normalizes each dataset independently."""
    payload = build_arity_hist_payload(
        pd.Series(["H2", "H2O"]),
        pd.Series(["LiFeO2", "NaCl"]),
        pd.Series(["SiO2", "Al2O3"]),
    )
    assert [series["label"] for series in payload["datasets"]] == [
        "MP",
        "MPtrj",
        "WBM",
    ]
    assert [series["color"] for series in payload["datasets"]] == list(SERIES_COLORS)
    assert all(
        sum(value for value in series["y"] if value is not None) == pytest.approx(1)
        for series in payload["datasets"]
    )


def test_build_mp_trj_hist_payload() -> None:
    """MPtrj payload combines scalar, nested, and structure-size histograms."""
    df_mp_trj = pd.DataFrame(
        {
            MbdKey.e_form_dft: [-1.0, -0.5],
            Key.forces: [
                [[1.0, -2.0, 0.0], [0.5, 0.0, 0.0]],
                [[3.0, 0.0, -1.0]],
            ],
            Key.stress_trace: [-0.1, 0.2],
            Key.magmoms: [[0.0, 1.0], [0.5]],
            Key.atom_nums: [[1, 8], [6, 1, 1]],
        }
    )
    payload = build_mp_trj_hist_payload(df_mp_trj)
    assert set(payload) == {"e-form", "forces", "stresses", "magmoms", "n-sites"}
    assert payload["n-sites"]["y"] == [0, 2]
    assert payload["n-sites"]["cumulative"] == [0.0, 1.0]


def test_build_route_and_comparison_element_counts() -> None:
    """Route counts and the MP/WBM comparison share occurrence inputs."""
    df_mp = pd.DataFrame({Key.formula: ["Li2O", "NaCl"]})
    df_wbm = pd.DataFrame(
        {Key.formula: ["LiF", "NaCl", "MgO", "AlN", "SiC"]},
        index=[f"wbm-{step}-1" for step in range(1, 6)],
    )
    df_mp_trj = pd.DataFrame({Key.formula: ["Li2O", "LiF"]})
    counts = build_route_element_counts(df_mp, df_wbm, df_mp_trj)

    assert {
        "mp-element-counts-by-occurrence",
        "wbm-element-counts-batch=1",
        "wbm-element-counts-arity=2",
        "mp-trj-element-counts-by-composition",
    } <= set(counts)

    payload = build_element_counts_payload(
        counts["mp-element-counts-by-occurrence"], df_wbm[Key.formula]
    )
    assert set(payload) == {"raw", "normalized"}
    assert all(
        {"label", "x", "y"} <= set(series)
        for variant in payload.values()
        for series in variant
    )


def test_available_static_builders_match_committed_payloads() -> None:
    """Builders reproducibly match committed payloads without the MPtrj cache."""
    from matbench_discovery.data import df_wbm

    df_mp = pd.read_csv(DataFiles.mp_energies.path, na_filter=False).set_index(
        Key.mat_id
    )
    route_mp_counts = pd.read_json(
        f"{ROOT}/site/src/routes/data/mp-element-counts-by-occurrence.json",
        typ="series",
    )
    expected_payloads = {
        "hist-wbm-e-form-per-atom": build_wbm_e_form_hist(df_wbm[MbdKey.e_form_wbm]),
        "hist-wbm-hull-dist": build_wbm_hull_dist_hist(df_wbm[MbdKey.each_true]),
        "mp-elemental-ref-energies": build_mp_elemental_ref_energies(
            mp_elem_ref_entries
        ),
        "spacegroup-sunbursts": build_spacegroup_sunbursts(
            df_mp[MbdKey.protostructure_spglib].str.split("_").str[2].astype(int),
            df_wbm[MbdKey.init_protostructure_spglib].str.split("_").str[2].astype(int),
        ),
        "element-counts-mp-vs-wbm": build_element_counts_payload(
            route_mp_counts,
            df_wbm.query(MbdKey.uniq_proto)[Key.formula],
        ),
    }
    for name, expected in expected_payloads.items():
        with gzip.open(f"{ROOT}/site/src/figs/{name}.json.gz", "rt") as file:
            assert json.load(file) == expected

    route_counts = build_route_element_counts(
        df_mp,
        df_wbm,
        pd.DataFrame({Key.formula: ["Li2O"]}),
    )
    for name, expected in route_counts.items():
        if name.startswith("mp-trj"):
            continue
        committed = pd.read_json(
            f"{ROOT}/site/src/routes/data/{name}.json", typ="series"
        )
        pd.testing.assert_series_equal(
            committed, expected, check_dtype=False, check_names=False
        )

    for count_mode in ("occurrence", "composition"):
        path = f"{ROOT}/site/src/routes/data/mp-trj-element-counts-by-{count_mode}.json"
        with open(path, encoding="utf-8") as file:
            payload = json.load(file)
        assert payload
        for value in payload.values():
            if value is not None:
                assert isinstance(value, int | float)
