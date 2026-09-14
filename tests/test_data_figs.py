"""Tests for model-independent data-page payload builders."""

import gzip
import json
from pathlib import Path

import ase.io
import h5py
import pandas as pd
import pymatviz as pmv
import pytest
from pymatgen.core import Composition
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


def test_export_benchmark_element_counts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Occurrence counts ignore atom multiplicity, frames, and duplicate functionals."""
    from scripts import export_data_fig_payloads as exporter

    monkeypatch.setattr(
        DataFiles, "path", property(lambda file: str(tmp_path / file.name))
    )
    monkeypatch.setattr(exporter, "ROUTE_DATA_DIR", str(tmp_path))
    ase.io.write(
        DataFiles.phonondb_pbe_103_structures.path,
        [ase.Atoms("H2O"), ase.Atoms("H2")],
        format="extxyz",
    )
    with h5py.File(DataFiles.dynamat_v1_0_md_trajectories.path, "w") as file:
        for name, numbers in {"water": [1, 1, 8], "hydrogen": [1, 1]}.items():
            group = file.create_group(name)
            group["atomic_numbers"] = numbers
            group.create_dataset("positions", shape=(10, len(numbers), 3))
            group.attrs["temperature_kelvin"] = 300
            group.attrs["dt_fs"] = 2
    with gzip.open(DataFiles.diatomics_dft_reference.path, "wt") as file:
        json.dump({"PBE": {"H-H": {}}, "r2SCAN": {"H-H": {}, "He-He": {}}}, file)

    exporter.export_benchmark_element_counts()
    assert json.loads((tmp_path / "benchmark-element-counts.json").read_text()) == {
        "phonondb": {"H": 2, "O": 1},
        "dynamat": {"H": 2, "O": 1},
        "diatomics": {"H": 1, "He": 1},
    }
    pd.DataFrame({"n_sites": [2, 2, 3]}).to_csv(DataFiles.wbm_summary.path, index=False)
    with gzip.open(DataFiles.phonondb_pbe_103_kappa_no_nac.path, "wt") as file:
        json.dump(
            [
                {"temperatures": [600, 300], "kappa_tot_avg": [100, 2]},
                {"temperatures": [300], "kappa_tot_avg": [8]},
            ],
            file,
        )
    exporter.export_benchmark_summary()
    summary = json.loads((tmp_path / "benchmark-reference-summary.json").read_text())
    assert summary["wbm_n_sites"] == {"x": [2, 3], "y": [2, 1]}
    assert summary["phonondb_kappa"] == [2, 8]
    assert summary["md_systems"] == [
        {
            "name": "hydrogen",
            "formula": "H2",
            "n_atoms": 2,
            "temperature": 300,
            "duration_ps": 0.018,
        },
        {
            "name": "water",
            "formula": "H2O",
            "n_atoms": 3,
            "temperature": 300,
            "duration_ps": 0.018,
        },
    ]


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


@pytest.mark.parametrize(
    "wbm_formulas",
    [
        ["LiF", "NaCl", "MgO", "AlN", "SiC"],
        ["Li2O", "Li2O", "Ca3(PO4)2", "H2", "NaN"],
        ["Fe0"],
        [],
    ],
)
def test_build_route_and_comparison_element_counts(wbm_formulas: list[str]) -> None:
    """Route counts and the MP/WBM comparison share occurrence inputs."""
    df_mp = pd.DataFrame({Key.formula: ["Li2O", "NaCl", "NaN"]})
    df_wbm = pd.DataFrame(
        {Key.formula: wbm_formulas},
        index=[f"wbm-{step}-1" for step in range(1, len(wbm_formulas) + 1)],
    )
    df_mp_trj = pd.DataFrame({Key.formula: ["Fe0.1Ni0.9O", "Fe0.1Ni0.9O", "LiF"]})
    counts = build_route_element_counts(df_mp, df_wbm, df_mp_trj)

    for dataset, formulas in (
        ("mp", df_mp[Key.formula].iloc[:2]),
        ("wbm", df_wbm[Key.formula]),
        ("mp-trj", df_mp_trj[Key.formula]),
    ):
        for count_mode in ("occurrence", "composition"):
            expected = pmv.count_elements(formulas, count_mode=count_mode)
            if dataset != "mp-trj":
                expected = expected.astype("Int64")
            pd.testing.assert_series_equal(
                counts[f"{dataset}-element-counts-by-{count_mode}"],
                expected,
                check_exact=True,
            )
    for batch, formula in enumerate(wbm_formulas, start=1):
        pd.testing.assert_series_equal(
            counts[f"wbm-element-counts-batch={batch}"],
            pmv.count_elements([formula]),
            check_exact=True,
        )
    for arity, formulas in df_wbm[Key.formula].groupby(
        df_wbm[Key.formula].map(Composition).map(len)
    ):
        pd.testing.assert_series_equal(
            counts[f"wbm-element-counts-arity={arity}"],
            pmv.count_elements(formulas),
            check_exact=True,
        )

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
