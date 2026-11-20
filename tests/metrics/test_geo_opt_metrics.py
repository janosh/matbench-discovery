"""Tests for geometry optimization metrics (RMSD AUC and symmetry changes vs DFT)."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from pymatviz.enums import Key

from matbench_discovery.data import make_file_ref, round_trip_yaml
from matbench_discovery.enums import MbdKey, Model
from matbench_discovery.metrics import geo_opt
from tests.utils import patch_model_metadata

_GEO = "models/alignn/alignn/2026-07-01-geo-opt-symprec=1e-2-moyo=0.4.2.csv.gz"
_GEO_NAN = "models/alignn/alignn/2026-07-02-geo-opt-symprec=1e-2-moyo=0.4.2.csv.gz"


@pytest.mark.parametrize(
    (
        "spg_diffs",
        "n_sym_ops_diffs",
        "expected_decrease",
        "expected_match",
        "expected_increase",
        "expected_lateral",
        "expected_mae",
    ),
    [
        # All matches
        ([0, 0, 0], [0, 0, 0], 0.0, 1.0, 0.0, 0.0, 0.0),
        # All decreases
        ([-1, -2, -1], [-2, -4, -2], 1.0, 0.0, 0.0, 0.0, 8 / 3),
        # All increases
        ([1, 2, 1], [2, 4, 2], 0.0, 0.0, 1.0, 0.0, 8 / 3),
        # Mixed cases
        ([0, -1, 1], [0, -2, 2], 1 / 3, 1 / 3, 1 / 3, 0.0, 4 / 3),
        # space group changed with the same number of symmetry ops (e.g. Pbca->Pnma)
        ([1, -1, 0], [0, 0, 0], 0, 1 / 3, 0, 2 / 3, 0.0),
        # Include some NaN values
        ([0, np.nan, 1], [0, np.nan, 2], 0.0, 0.5, 0.5, 0.0, 1.0),
        # NaN in the middle: 5 valid rows with 1 match, 1 decrease, 3 increases
        (
            [0, 1, np.nan, -1, 2, 3],
            [0, 2, np.nan, -2, 4, 6],
            1 / 5,
            1 / 5,
            3 / 5,
            0.0,
            14 / 5,
        ),
        # no valid symmetry rows: every symmetry metric is NaN, not a 0/0 crash
        ([np.nan] * 2, [np.nan] * 2, np.nan, np.nan, np.nan, np.nan, np.nan),
    ],
)
def test_calc_geo_opt_metrics_parametrized(
    spg_diffs: list[float],
    n_sym_ops_diffs: list[float],
    expected_decrease: float,
    expected_match: float,
    expected_increase: float,
    expected_lateral: float,
    expected_mae: float,
) -> None:
    """Test calc_geo_opt_metrics with various symmetry difference patterns."""
    # RMSD NaN pattern deliberately differs from the spg NaN pattern: n_structures
    # must count valid symmetry rows only, while the unmatched (NaN) RMSD and the
    # reference ID the model lacks both score 0 in rmsd_auc
    ids = [f"wbm-{idx}" for idx in range(len(spg_diffs))]
    rmsds = [0.01] * (len(spg_diffs) - 1) + [np.nan]
    df_geo_opt = pd.DataFrame(
        {
            MbdKey.structure_rmsd_vs_dft: rmsds,
            MbdKey.spg_num_diff: spg_diffs,
            MbdKey.n_sym_ops_diff: n_sym_ops_diffs,
        },
        index=ids,
    )
    # an ID outside the test set is dropped, one the model lacks counts as a failure
    df_geo_opt.loc["not-in-test-set"] = [0.0, 0, 0]
    ref_ids = pd.Index([*ids, "missing-prediction"])

    results = geo_opt.calc_geo_opt_metrics(df_geo_opt, ref_ids)

    for key, expected in (
        (Key.symmetry_decrease, expected_decrease),
        (Key.symmetry_match, expected_match),
        (Key.symmetry_increase, expected_increase),
        (MbdKey.symmetry_lateral, expected_lateral),
        (Key.n_sym_ops_mae, expected_mae),
    ):
        assert results[str(key)] == pytest.approx(expected, nan_ok=True), key
    buckets = (Key.symmetry_decrease, Key.symmetry_match, Key.symmetry_increase)
    if results[str(Key.n_structures)]:  # the 4 buckets partition valid structures
        fractions = [results[str(key)] for key in (*buckets, MbdKey.symmetry_lateral)]
        assert sum(fractions) == pytest.approx(1)
    # n_structures should be the number of non-NaN spg_diff values
    assert results[str(Key.n_structures)] == np.count_nonzero(pd.notna(spg_diffs))
    expected_auc = (len(spg_diffs) - 1) * (1 - 0.01 / 0.05) / len(ref_ids)
    assert results[str(MbdKey.rmsd_auc)] == pytest.approx(expected_auc, rel=1e-12)


@pytest.mark.parametrize(
    ("rmsds", "ref_ids", "expected"),
    [
        ([0.0, 0.0], ["a", "b"], 1.0),  # perfect reproduction of DFT
        ([0.025, 0.05], ["a", "b"], 0.25),  # half credit, and 0 at the cutoff
        ([0.2, np.nan], ["a", "b"], 0.0),  # beyond the cutoff and unmatched score 0
        ([0.0], ["a", "b"], 0.5),  # the missing reference ID scores 0
        ([0.0, 0.01], ["a"], 1.0),  # predictions outside the test set are ignored
    ],
)
def test_rmsd_auc(rmsds: list[float], ref_ids: list[str], expected: float) -> None:
    """rmsd_auc averages max(1 - rmsd/cutoff, 0) over all reference IDs."""
    series = pd.Series(rmsds, index=list("ab")[: len(rmsds)])
    assert geo_opt.rmsd_auc(series, pd.Index(ref_ids)) == pytest.approx(
        expected, rel=1e-12
    )


def test_rmsd_auc_equals_cdf_integral() -> None:
    """The closed form equals the normalized integral of the empirical CDF."""
    rng = np.random.default_rng(seed=0)
    rmsds = pd.Series(rng.exponential(0.02, 500))
    rmsds[rng.choice(500, 40, replace=False)] = np.nan
    ref_ids = pd.RangeIndex(550)  # 50 missing predictions
    grid = np.linspace(0, geo_opt.RMSD_AUC_CUTOFF, 2_000_001)
    cdf = np.searchsorted(np.sort(rmsds.dropna()), grid, side="right") / len(ref_ids)
    numeric = np.trapezoid(cdf, grid) / geo_opt.RMSD_AUC_CUTOFF
    # each CDF jump (total height <= 1) costs the trapezoid at most dx / 2 of area,
    # so the normalized error is <= dx / (2 * cutoff) = 2.5e-7 on this grid
    assert geo_opt.rmsd_auc(rmsds, ref_ids) == pytest.approx(numeric, abs=2.5e-7)


@pytest.mark.parametrize(
    ("rmsds", "ref_ids", "match"),
    [
        (pd.Series([0.0]), [], "at least one reference"),
        # a CSV read without index_col: integer index never overlaps material IDs
        (pd.Series([0.0, 0.0]), ["wbm-1", "wbm-2"], "no RMSD index entries"),
        (pd.Series([0.0], index=["wbm-1"]), ["wbm-1", "wbm-1"], "duplicate ref_ids"),
        (pd.Series([0.0, 0.1], index=["wbm-1"] * 2), ["wbm-1"], "duplicate rmsds"),
        (pd.Series([True], index=["wbm-1"]), ["wbm-1"], "must be floats"),
        (pd.Series([-0.05], index=["wbm-1"]), ["wbm-1"], "negative RMSDs"),
    ],
)
def test_rmsd_auc_rejects_bad_input(
    rmsds: pd.Series, ref_ids: list[str], match: str
) -> None:
    """Malformed inputs raise instead of returning a plausible-looking score."""
    with pytest.raises(ValueError, match=match):
        geo_opt.rmsd_auc(rmsds, pd.Index(ref_ids, dtype=object))


_GEO_REF_WITH_META = make_file_ref(
    _GEO, url="https://figshare.com/files/1", size=7, md5="a" * 32
)
_FULL_METRICS = {
    MbdKey.rmsd_auc: 0.8,
    Key.n_sym_ops_mae: 0.2,
    Key.symmetry_decrease: 0.3,
    Key.symmetry_match: 0.4,
    Key.symmetry_increase: 0.0,
    MbdKey.symmetry_lateral: 0.0,
    Key.n_structures: 0,
}
_FULL_BLOCK = {
    Key.n_sym_ops_mae: 0.2,
    Key.symmetry_decrease: 0.3,
    Key.symmetry_match: 0.4,
    Key.symmetry_increase: 0.0,
    MbdKey.symmetry_lateral: 0.0,
    Key.n_structures: 0,
}


@pytest.mark.parametrize(
    ("metrics_data", "expected_block", "analysis_file_path", "existing_ref"),
    [
        (
            _FULL_METRICS,
            {**_FULL_BLOCK, "analysis_file": make_file_ref(_GEO)},
            _GEO,
            None,
        ),
        (
            {
                MbdKey.rmsd_auc: 0.0,
                Key.n_sym_ops_mae: float("nan"),
                Key.symmetry_decrease: 0.0,
                Key.symmetry_match: 0.0,
                Key.symmetry_increase: 0.0,
                MbdKey.symmetry_lateral: 0.0,
                Key.n_structures: 0,
            },
            {
                Key.n_sym_ops_mae: float("nan"),
                Key.symmetry_decrease: 0.0,
                Key.symmetry_match: 0.0,
                Key.symmetry_increase: 0.0,
                MbdKey.symmetry_lateral: 0.0,
                Key.n_structures: 0,
                # renamed analysis file drops the stale url/size/md5
                "analysis_file": make_file_ref(_GEO_NAN),
            },
            _GEO_NAN,
            _GEO_REF_WITH_META,
        ),
        # same file name (metrics re-derived from a cached CSV) keeps its Figshare ref
        (
            _FULL_METRICS,
            {**_FULL_BLOCK, "analysis_file": _GEO_REF_WITH_META},
            _GEO,
            _GEO_REF_WITH_META,
        ),
    ],
    ids=["fresh", "nan-renamed-file", "cached-keeps-ref"],
)
def test_write_geo_opt_metrics_to_yaml(
    metrics_data: dict[MbdKey | Key, float],
    expected_block: dict[MbdKey | Key | str, float | str | None],
    analysis_file_path: str,
    existing_ref: dict[str, object] | None,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Test saving geometry optimization metrics to YAML files with edge cases."""
    symprec = 1e-2
    symprec_key = "symprec=1e-2"

    real_yaml_path = Path(Model["alignn"].yaml_path)
    original = real_yaml_path.read_bytes()
    yaml_path = tmp_path / "model.yml"
    prior_block = {"analysis_file": existing_ref} if existing_ref else {}
    with yaml_path.open("w", encoding="utf-8") as file:
        round_trip_yaml.dump({"metrics": {"geo_opt": {symprec_key: prior_block}}}, file)
    patch_model_metadata(monkeypatch, Model["alignn"], yaml_path=str(yaml_path))

    geo_opt.write_metrics_to_yaml(
        pd.DataFrame([metrics_data]), Model["alignn"], symprec, analysis_file_path
    )
    with yaml_path.open(encoding="utf-8") as file:
        actual_block = round_trip_yaml.load(file)["metrics"]["geo_opt"][symprec_key]
    assert set(actual_block) == set(expected_block)
    for key, expected_val in expected_block.items():
        value = actual_block[key]
        if key == "analysis_file":
            assert value == expected_val
        else:
            assert value == pytest.approx(expected_val, nan_ok=True)
    # RMSD is symprec-invariant: its score sits once at task level
    geo_opt_block = round_trip_yaml.load(yaml_path.read_text())["metrics"]["geo_opt"]
    assert geo_opt_block[str(MbdKey.rmsd_auc)] == metrics_data[MbdKey.rmsd_auc]
    assert "rmsd" not in actual_block
    assert real_yaml_path.read_bytes() == original
