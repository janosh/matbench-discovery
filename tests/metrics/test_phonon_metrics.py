"""Tests for thermal conductivity metrics."""

from collections.abc import Mapping
from functools import partial
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose
from numpy.typing import NDArray
from pymatviz.enums import Key

from matbench_discovery import data as mbd_data
from matbench_discovery.data import make_file_ref
from matbench_discovery.enums import MbdKey, Model
from matbench_discovery.metrics import phonons as phonon_metrics

_KAPPA_MODEL_DIR = "models/mace/mace-mp-0"
KAPPA_PRED = f"{_KAPPA_MODEL_DIR}/2026-07-01-phonons-kappa-103.json.gz"
KAPPA_FORCE = f"{_KAPPA_MODEL_DIR}/2026-07-01-phonons-kappa-103-forces.json.gz"
KAPPA_RUN_INFO = f"{_KAPPA_MODEL_DIR}/2026-07-01-phonons-kappa-103-run-info.json"
KAPPA_PRED_EXISTING = f"{_KAPPA_MODEL_DIR}/2025-01-01-phonons-kappa-103.json.gz"
KAPPA_FORCE_NEW = f"{_KAPPA_MODEL_DIR}/2026-07-02-phonons-kappa-103-forces.json.gz"
_GEO_OPT_FILE = "models/alignn/alignn/2026-07-01-geo-opt-symprec=1e-2-moyo=0.4.2.csv.gz"
KAPPA_METRICS: dict[str, float | None] = {
    "srme": 0.25,
    "sre": 0.125,
    "srd": -0.05,
    "failure_rate": 0.1,
    "imaginary_mode_rate": 0.05,
    "spectrum_w1": 0.2,
}


@pytest.fixture
def df_pred() -> pd.DataFrame:
    """Mock DataFrame with ML predictions."""
    return pd.DataFrame(
        {
            MbdKey.kappa_tot_rta: [np.diag([1, 2, 3]), 2 * np.eye(3)],
            MbdKey.kappa_tot_avg: [2, 2],  # average of diagonal elements
            MbdKey.mode_kappa_tot_rta: [np.diag([0.5, 1, 1.5]), np.eye(3)],
            Key.mode_weights: [np.array([1]), np.array([1])],
            Key.has_imag_ph_modes: [False, False],
            Key.final_spg_num: [1, 1],
            Key.init_spg_num: [1, 1],
        }
    )


@pytest.fixture
def df_true() -> pd.DataFrame:
    """Mock DataFrame with DFT reference values."""
    return pd.DataFrame(
        {
            MbdKey.kappa_tot_rta: [np.diag([1, 2, 3]), np.diag([2, 2, 2])],
            MbdKey.kappa_tot_avg: [2, 2],
            MbdKey.mode_kappa_tot_rta: [np.diag([0.5, 1, 1.5]), np.eye(3)],
            Key.mode_weights: [np.array([1]), np.array([1])],
        }
    )


@pytest.fixture
def df_minimal() -> pd.DataFrame:
    """Minimal DataFrame with required columns."""
    return pd.DataFrame(
        {
            MbdKey.kappa_tot_rta: [np.eye(3)],
            MbdKey.kappa_tot_avg: [1.0],
            Key.mode_weights: [np.array([1])],
            MbdKey.mode_kappa_tot_rta: [np.eye(3)],
        }
    )


@pytest.fixture
def series_multi_temp() -> pd.Series:
    """Mock Series with multi-temperature data."""
    return pd.Series(
        {
            MbdKey.kappa_tot_avg: np.array([2.0, 1.5, 1.0]),
            MbdKey.kappa_tot_rta: np.array([2 * np.eye(3), 1.5 * np.eye(3), np.eye(3)]),
            MbdKey.mode_kappa_tot_rta: np.array(
                [[2 * np.eye(3)], [1.5 * np.eye(3)], [np.eye(3)]]
            ),
            Key.mode_weights: np.ones(1),
        }
    )


@pytest.mark.parametrize(
    "tensor,expected,description",
    [
        (np.diag([1, 2, 3]), 2, "diagonal tensor"),
        (
            [[1, 0.1, 0], [0.1, 2, 0], [0, 0, 3]],
            2,
            "non-diagonal tensor",
        ),
        (np.zeros((3, 3)), 0, "zero tensor"),
        (-np.diag([1, 2, 3]), -2, "negative tensor"),
        # Voigt 6-vectors [xx, yy, zz, yz, xz, xy]: the format kappa_tot_rta is
        # stored in by calc_kappa.py (regression test: a rewrite once dropped Voigt
        # support, silently corrupting CHGNet/M3GNet kappa_SRME to 2.0 on re-eval)
        (np.array([1.0, 2, 3, 0.5, 0.5, 0.5]), 2, "Voigt 6-vector"),
        (
            np.array([[1.0, 2, 3, 0, 0, 0], [2, 4, 6, 0, 0, 0]]),
            [2, 4],
            "multi-temperature Voigt",
        ),
        (
            np.array([np.eye(3), 2 * np.eye(3)]),
            [1, 2],
            "multi-temperature tensors",
        ),
        (np.diag([np.nan, 2, 3]), np.nan, "NaN tensor"),
        (np.asarray(None, dtype=float), np.nan, "failed scalar"),
    ],
)
def test_calculate_kappa_avg(
    tensor: object, expected: object, description: str
) -> None:
    """Test calculation of average thermal conductivity with various inputs."""
    avg = phonon_metrics.calculate_kappa_avg(np.asarray(tensor))
    assert_allclose(avg, np.asarray(expected), rtol=1e-6, err_msg=description)


@pytest.mark.parametrize("n_levels", [5, 11, 101])
def test_weighted_quantiles_uniform_matches_numpy_inverted_cdf(n_levels: int) -> None:
    """Uniform-weight quantiles equal numpy's inverted_cdf method exactly."""
    np_rng = np.random.default_rng(seed=0)
    values = np_rng.normal(size=500)
    levels = np.linspace(0, 1, n_levels)
    result = phonon_metrics.weighted_quantiles(values, None, levels)
    expected = np.quantile(values, levels, method="inverted_cdf")
    assert_allclose(result, expected, rtol=0, atol=0)  # bit-identical


def test_weighted_quantiles_integer_weights_equal_repetition() -> None:
    """Integer weights are equivalent to repeating samples that many times."""
    values = np.array([0.0, 1.0, 2.0, 4.0, 8.0, 16.0])
    weights = np.array([0.0, 1.0, 3.0, 2.0, 1.0, 0.0])
    repeated = np.repeat(values, weights.astype(int))
    levels = np.linspace(0, 1, 29)
    result = phonon_metrics.weighted_quantiles(values, weights, levels)
    expected = phonon_metrics.weighted_quantiles(repeated, None, levels)
    assert_allclose(result, expected, rtol=0, atol=0)  # bit-identical


def test_weighted_quantiles_drops_nans_and_rejects_empty() -> None:
    """NaN values are dropped with their weights; all-NaN input raises."""
    values = np.array([1.0, np.nan, 3.0])
    weights = np.array([1.0, 100.0, 1.0])
    levels = np.array([0.0, 0.75, 1.0])
    result = phonon_metrics.weighted_quantiles(values, weights, levels)
    # NaN (and its weight 100) dropped -> remaining CDF is [0.5, 1.0]
    assert_allclose(result, [1.0, 3.0, 3.0], rtol=0, atol=0)
    with pytest.raises(ValueError, match="no finite values"):
        phonon_metrics.weighted_quantiles(np.array([np.nan]), None, levels)
    with pytest.raises(ValueError, match="!="):
        phonon_metrics.weighted_quantiles(values, np.ones(2), levels)


@pytest.mark.parametrize(
    "weights,levels,match",
    [
        (np.array([0.0, 0.0, 0.0]), np.array([0.5]), "positive"),
        (np.array([1.0, -1.0, 1.0]), np.array([0.5]), "non-negative"),
        (np.array([1.0, np.nan, 1.0]), np.array([0.5]), "finite"),
        (np.array([1.0, 1.0, 1.0]), np.array([-0.1]), r"\[0, 1\]"),
        (np.array([1.0, 1.0, 1.0]), np.array([1.1]), r"\[0, 1\]"),
        (np.array([1.0, 1.0, 1.0]), np.array([np.nan]), r"\[0, 1\]"),
    ],
)
def test_weighted_quantiles_rejects_invalid_weights_and_levels(
    weights: NDArray[np.float64], levels: NDArray[np.float64], match: str
) -> None:
    """Invalid weights and quantile levels raise instead of producing endpoints."""
    with pytest.raises(ValueError, match=match):
        phonon_metrics.weighted_quantiles(np.array([1.0, 2.0, 3.0]), weights, levels)


@pytest.mark.parametrize(
    ("weights", "expected_len"),
    [([1.0, 2.0, 1.0, 1.0], 12), (None, None), ([1.0, 2.0], None)],
)
def test_mode_weights_for_freqs(
    weights: list[float] | None, expected_len: int | None
) -> None:
    """q-point weights repeat per mode iff their length matches the q-point axis."""
    n_modes = 3
    ph_freqs = np.ones((4, n_modes))
    q_weights = None if weights is None else np.array(weights)
    result = phonon_metrics.mode_weights_for_freqs(ph_freqs, q_weights)
    if expected_len is None:
        assert result is None
    else:
        assert result is not None
        assert weights is not None
        assert len(result) == expected_len
        assert_allclose(result[:n_modes], weights[0], rtol=0, atol=0)


@pytest.mark.parametrize(
    ("overrides", "series_srme", "censored"),
    [
        # imaginary modes and broken symmetry are censored at the DataFrame level,
        # so the per-material series SRME of the otherwise exact prediction is 0
        ({Key.has_imag_ph_modes: True}, 0, True),
        ({Key.final_spg_num: 2}, 0, True),
        ({MbdKey.kappa_tot_avg: np.nan}, 2, True),
        ({MbdKey.kappa_tot_avg: np.inf}, 2, True),
        ({MbdKey.kappa_tot_avg: -1}, 2, True),
        # zero predicted kappa is a valid exact-2 score, not a censored failure
        (
            {MbdKey.kappa_tot_avg: 0, MbdKey.mode_kappa_tot_rta: np.zeros((3, 3))},
            2,
            False,
        ),
    ],
)
def test_calc_kappa_srme_invalid_predictions_score_two(
    df_pred: pd.DataFrame,
    df_true: pd.DataFrame,
    overrides: dict[str, object],
    series_srme: float,
    censored: bool,
) -> None:
    """Imaginary modes, broken symmetry, invalid and zero kappa get maximum SRME."""
    modified = df_pred.copy()
    for column, value in overrides.items():
        modified[column] = modified[column].astype(object)
        modified.at[0, column] = value  # noqa: PD008 (.loc broadcasts arrays)
    # match the reference total to its mode-kappa average (diag 0.5, 1, 1.5) so an
    # all-zero prediction scores exactly 2; SRME numerators only use mode kappas
    reference = df_true.copy()
    reference.loc[0, MbdKey.kappa_tot_avg] = 1
    srme_values, censored_flags = phonon_metrics.calc_kappa_srme_dataframes(
        modified, reference
    )
    assert (srme_values[0], censored_flags[0]) == (2, censored)
    series_srmes = phonon_metrics.calc_kappa_srme(modified.loc[0], reference.loc[0])
    assert list(series_srmes) == [series_srme]


def test_calc_kappa_srme_single_material() -> None:
    """Test SRME calculation for a single material."""
    ml_data = pd.Series(
        {
            MbdKey.kappa_tot_avg: np.array([2.0]),
            MbdKey.kappa_tot_rta: 2 * np.eye(3),
            MbdKey.mode_kappa_tot_rta: np.eye(3),
            Key.mode_weights: np.array([1.0]),
        }
    )
    dft_data = ml_data.copy()

    assert phonon_metrics.calc_kappa_srme(ml_data, dft_data)[0] == 0
    ml_data[MbdKey.mode_kappa_tot_avg] = dft_data[MbdKey.mode_kappa_tot_avg] = np.nan
    assert phonon_metrics.calc_kappa_srme(ml_data, dft_data)[0] == 0

    ml_data[MbdKey.kappa_tot_avg] = np.array([3.0])
    ml_data[MbdKey.mode_kappa_tot_rta] = [1.5 * np.eye(3)]
    assert phonon_metrics.calc_kappa_srme(ml_data, dft_data)[0] > 0


def test_calc_kappa_metrics_with_different_values(
    df_pred: pd.DataFrame, df_true: pd.DataFrame
) -> None:
    """Per-material SRD/SRE/SRME compare predicted tensors with reference scalars."""
    df_pred_copy = df_pred.copy()
    df_pred_copy[MbdKey.kappa_tot_rta] = [2 * np.diag([1, 2, 3]), 4 * np.eye(3)]
    df_pred_copy[MbdKey.mode_kappa_tot_rta] = [2 * np.diag([1, 2, 3]), 4 * np.eye(3)]

    df_out = phonon_metrics.calc_kappa_metrics_from_dfs(df_pred_copy, df_true)
    pd.testing.assert_index_equal(df_out.index, df_pred.index)
    assert set(df_out) == {
        Key.sre,
        Key.srme,
        Key.final_spg_num,
        Key.has_imag_ph_modes,
        Key.init_spg_num,
        MbdKey.kappa_tot_avg,
        MbdKey.kappa_tot_rta,
        MbdKey.mode_kappa_tot_rta,
        Key.mode_weights,
        Key.srd,
        phonon_metrics.SRME_CENSORED,
        MbdKey.true_kappa_tot_avg,
    }
    # kappa_tot_avg is recomputed from the tensors: 4 vs reference 2 for both rows
    assert_allclose(df_out[Key.srd], [2 / 3, 2 / 3], rtol=1e-12)
    assert_allclose(df_out[Key.sre], [2 / 3, 2 / 3], rtol=1e-12)
    assert df_out[Key.srme].mean() == pytest.approx(1)
    assert df_out[phonon_metrics.SRME_CENSORED].tolist() == [False, False]


def test_calc_kappa_srme_temperature_dependence(series_multi_temp: pd.Series) -> None:
    """Test SRME calculation with temperature-dependent conductivities."""
    ml_data = series_multi_temp.copy()
    dft_data = series_multi_temp.copy()
    assert list(phonon_metrics.calc_kappa_srme(ml_data, dft_data)) == [0, 0, 0]
    # Make DFT values half of ML predictions without mutating shared arrays.
    dft_data[MbdKey.kappa_tot_avg] = ml_data[MbdKey.kappa_tot_avg] / 2
    dft_data[MbdKey.kappa_tot_rta] = ml_data[MbdKey.kappa_tot_rta] / 2
    dft_data[MbdKey.mode_kappa_tot_rta] = ml_data[MbdKey.mode_kappa_tot_rta] / 2

    kappa_srmes = phonon_metrics.calc_kappa_srme(ml_data, dft_data)
    assert list(kappa_srmes) == pytest.approx([2 / 3, 2 / 3, 2 / 3])

    # Dataframe metrics reject ambiguous multi-temperature rows.
    ml_data[Key.has_imag_ph_modes] = False
    ml_data[Key.final_spg_num] = 1
    ml_data[Key.init_spg_num] = 1

    with pytest.raises(ValueError, match="Reference kappa"):
        phonon_metrics.calc_kappa_srme_dataframes(
            pd.DataFrame([ml_data]), pd.DataFrame([dft_data])
        )


@pytest.mark.parametrize("invalid_value", [np.nan, np.inf, -1.0])
def test_calc_kappa_srme_rejects_invalid_reference(invalid_value: float) -> None:
    """Invalid reference totals abort scoring instead of penalizing a model."""
    prediction = pd.Series(
        {
            MbdKey.kappa_tot_avg: np.array([1.0]),
            MbdKey.mode_kappa_tot_rta: np.eye(3),
            Key.mode_weights: np.ones(3),
        }
    )
    reference = prediction.copy()
    reference[MbdKey.kappa_tot_avg] = np.array([invalid_value])
    with pytest.raises(ValueError, match="Reference kappa totals"):
        phonon_metrics.calc_kappa_srme(prediction, reference)
    if np.isnan(invalid_value):
        reference[MbdKey.kappa_tot_avg] = np.array([1.0])
        reference[MbdKey.mode_kappa_tot_rta] = np.full((3, 3), np.nan)
        prediction[Key.has_imag_ph_modes] = True
        with pytest.raises(ValueError, match="Reference mode"):
            phonon_metrics.calc_kappa_srme_dataframes(
                pd.DataFrame([prediction]), pd.DataFrame([reference])
            )
        prediction[MbdKey.kappa_tot_avg] = np.array([np.nan])
        with pytest.raises(ValueError, match="Reference mode"):
            phonon_metrics.calc_kappa_srme(prediction, reference)


def test_calc_kappa_metrics_from_dfs_symmetry(df_minimal: pd.DataFrame) -> None:
    """Distinguish valid exact-2 SRME from imaginary-mode and symmetry failures."""
    df_pred = pd.concat([df_minimal] * 3, ignore_index=True).copy()
    df_pred[MbdKey.kappa_tot_rta] = [np.zeros((3, 3)), np.eye(3), np.eye(3)]
    df_pred[MbdKey.mode_kappa_tot_rta] = [np.zeros((3, 3)), np.eye(3), np.eye(3)]
    df_pred[Key.has_imag_ph_modes] = [False, True, False]
    df_pred[Key.final_spg_num] = [1, 1, 2]
    df_pred[Key.init_spg_num] = [1, 1, 1]

    df_true = pd.concat([df_minimal] * 3, ignore_index=True).copy()
    df_true[Key.init_spg_num] = [1, 1, 1]

    result = phonon_metrics.calc_kappa_metrics_from_dfs(df_pred, df_true)
    assert result[Key.srme].eq(2).tolist() == [True, True, True]
    assert result[phonon_metrics.SRME_CENSORED].tolist() == [False, True, True]


@pytest.mark.parametrize(
    ("reference_spg_col", "expected_failures"),
    [
        (Key.init_spg_num, [False, True]),
        (None, [False, False]),
    ],
)
def test_calc_kappa_srme_dataframes_missing_init_spg_uses_reference(
    df_minimal: pd.DataFrame,
    reference_spg_col: str | None,
    expected_failures: list[bool],
) -> None:
    """Missing prediction symmetry uses canonical or absent references."""
    df_pred = pd.concat([df_minimal] * 2, ignore_index=True).copy()
    df_pred[Key.final_spg_num] = [225, 186]
    df_true = pd.concat([df_minimal] * 2, ignore_index=True).copy()
    if reference_spg_col is not None:
        df_true[reference_spg_col] = [225, 225]

    srme_values, _censored = phonon_metrics.calc_kappa_srme_dataframes(df_pred, df_true)
    assert [value == 2 for value in srme_values] == expected_failures


def test_calc_kappa_srme_dataframes_penalizes_float_typed_imaginary_flag(
    df_minimal: pd.DataFrame,
) -> None:
    """Imaginary-mode rows score 2 even when the flag deserializes as 1.0/NaN.

    A has_imag_ph_modes column with missing values becomes float dtype, so an
    `is True` identity check would silently skip the penalty for those files.
    """
    df_pred = pd.concat([df_minimal] * 3, ignore_index=True).copy()
    df_pred[Key.has_imag_ph_modes] = [1.0, 0.0, np.nan]
    df_true = pd.concat([df_minimal] * 3, ignore_index=True).copy()

    srme_values, _censored = phonon_metrics.calc_kappa_srme_dataframes(df_pred, df_true)
    assert [value == 2 for value in srme_values] == [True, False, False]


def test_evaluate_kappa_predictions_returns_all_aggregate_metrics(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Aggregate scoring includes signed bias, failures, imaginary modes, and W1."""
    from matbench_discovery import phonons

    material_ids = ["mp-1", "mp-2", "mp-3"]
    failed_tensor = np.full((3, 3), np.nan)
    reference = pd.DataFrame(
        {
            MbdKey.kappa_tot_rta: [np.eye(3)] * 3,
            MbdKey.kappa_tot_avg: [1.0] * 3,
            MbdKey.mode_kappa_tot_rta: [np.eye(3)] * 3,
            Key.mode_weights: [np.ones(1)] * 3,
            Key.ph_freqs: [np.array([[0.0, 1.0, 2.0]])] * 3,
        },
        index=material_ids,
    )
    predictions = pd.DataFrame(
        {
            MbdKey.kappa_tot_rta: [
                np.eye(3),
                3 * np.eye(3),
                failed_tensor,
            ],
            MbdKey.mode_kappa_tot_rta: [
                np.eye(3),
                3 * np.eye(3),
                failed_tensor,
            ],
            Key.mode_weights: [np.ones(1)] * 3,
            Key.has_imag_ph_modes: [False, np.nan, True],
            Key.final_spg_num: [1, 1, 1],
            Key.init_spg_num: [1, 1, 1],
            Key.ph_freqs: [
                np.array([[0.0, 1.0, 2.0]]),
                np.array([[1.0, 2.0, 3.0]]),
                np.array([[np.nan, 1.0, 2.0]]),
            ],
        },
        index=material_ids,
    )

    def read_reference(_path: str) -> pd.DataFrame:
        """Return the fixed three-material reference frame."""
        return reference

    monkeypatch.setattr(phonons, "read_kappa_json", read_reference)
    metrics = phonon_metrics.evaluate_kappa_predictions(predictions)

    assert metrics == pytest.approx(
        {
            "srme": 1.0,
            "sre": 1.0,
            "srd": -1 / 3,
            "failure_rate": 1 / 3,
            "imaginary_mode_rate": 1 / 3,
            "spectrum_w1": 0.5,
        }
    )
    # an absent imaginary-mode column must fail loud, not score a 0% rate
    with pytest.raises(KeyError, match="has_imag_ph_modes"):
        phonon_metrics.evaluate_kappa_predictions(
            predictions.drop(columns=[Key.has_imag_ph_modes])
        )
    predictions[Key.ph_freqs] = [None] * len(predictions)
    no_spectrum_metrics = phonon_metrics.evaluate_kappa_predictions(predictions)
    assert no_spectrum_metrics["spectrum_w1"] is None
    for key in (MbdKey.kappa_tot_rta, MbdKey.mode_kappa_tot_rta):
        predictions[key] = [failed_tensor] * len(predictions)
    predictions[Key.has_imag_ph_modes] = True
    failed_metrics = phonon_metrics.evaluate_kappa_predictions(predictions)
    assert failed_metrics["srd"] == -2
    assert failed_metrics["sre"] == 2
    assert failed_metrics["failure_rate"] == 1
    assert failed_metrics["imaginary_mode_rate"] == 1


@pytest.mark.parametrize(
    "prediction_freqs",
    [None, [], [np.nan], [[np.nan, 1.0]]],
)
def test_phonon_spectrum_w1_rejects_unusable_predictions(
    prediction_freqs: object,
) -> None:
    """Missing, malformed, or non-finite predicted spectra have no W1 value."""
    prediction = pd.Series({Key.ph_freqs: prediction_freqs})
    reference = pd.Series({Key.ph_freqs: [[0.0, 1.0]]})
    assert phonon_metrics.phonon_spectrum_w1(prediction, reference) is None


@pytest.mark.parametrize("q_weights", [["bad"], [np.nan], [-1.0], [0.0]])
def test_phonon_spectrum_w1_falls_back_from_invalid_weights(
    q_weights: object,
) -> None:
    """Malformed q-point weights fall back to a uniform spectrum comparison."""
    prediction = pd.Series({Key.ph_freqs: [[0.0, 1.0]], "weights": q_weights})
    reference = pd.Series({Key.ph_freqs: [[0.0, 1.0]]})
    assert phonon_metrics.phonon_spectrum_w1(prediction, reference) == 0


def test_kappa_analysis_uses_canonical_spacegroups() -> None:
    """Analysis reads normalized initial space-group numbers."""
    from scripts.model_figs.kappa_103_analysis import canonical_spacegroup_numbers

    df_reference = pd.DataFrame({Key.init_spg_num: [225]}, index=["mp-1"])
    assert canonical_spacegroup_numbers(df_reference, ["mp-1"]) == [225]


@pytest.mark.parametrize("spacegroup_num", [None, np.nan, 0, 231, 1.5, True, "bad"])
def test_kappa_analysis_rejects_invalid_spacegroups(spacegroup_num: object) -> None:
    """Analysis rejects space groups that cannot drive crystal-system plots."""
    from scripts.model_figs.kappa_103_analysis import canonical_spacegroup_numbers

    df_reference = pd.DataFrame({Key.init_spg_num: [spacegroup_num]}, index=["mp-1"])
    with pytest.raises(ValueError, match=r"Invalid canonical space-group.*mp-1"):
        canonical_spacegroup_numbers(df_reference, ["mp-1"])


def test_kappa_analysis_uses_shared_spectrum_helpers(
    df_pred: pd.DataFrame, df_true: pd.DataFrame
) -> None:
    """The analysis payload reports the canonical metric's spectrum W1 per material."""
    from scripts.model_figs.kappa_103_analysis import model_payload

    material_ids = ["mp-1", "mp-2"]
    df_pred = df_pred.set_axis(material_ids)
    df_true = df_true.set_axis(material_ids)
    df_pred[Key.ph_freqs] = [np.array([[0.0, 1.0, 2.0]]), np.array([[0.5, 1.0, 2.0]])]
    df_true[Key.ph_freqs] = [np.array([[0.0, 1.0, 2.0]])] * len(material_ids)
    payload = model_payload(
        Model["mace_mp_0"], __file__, df_pred, df_true, material_ids
    )

    expected_w1 = [
        phonon_metrics.phonon_spectrum_w1(df_pred.loc[mat_id], df_true.loc[mat_id])
        for mat_id in material_ids
    ]
    assert expected_w1 == pytest.approx([0, 0.5 / 3])  # one of 3 modes shifted 0.5
    assert payload["srme_censored"] == [False, False]
    assert payload["freq_w1"] == pytest.approx(expected_w1, abs=5e-5)
    assert payload["freq_w1_mean"] == pytest.approx(np.mean(expected_w1), abs=1e-4)


def test_evaluate_kappa_predictions_rejects_duplicate_ids(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Duplicate prediction IDs cannot silently receive extra metric weight."""
    from matbench_discovery import phonons

    def read_reference(_path: str) -> pd.DataFrame:
        """Return a two-material reference frame."""
        return pd.DataFrame(index=["mp-1", "mp-2"])

    monkeypatch.setattr(phonons, "read_kappa_json", read_reference)
    with pytest.raises(ValueError, match="duplicate material IDs"):
        phonon_metrics.evaluate_kappa_predictions(pd.DataFrame(index=["mp-1", "mp-1"]))


def kappa_yaml_model(tmp_path: Path, kappa_103: Mapping[str, object]) -> Model:
    """Model stub whose temp YAML holds the given metrics.phonons.kappa_103."""
    yaml_path = f"{tmp_path}/model.yml"
    with open(yaml_path, mode="w", encoding="utf-8") as file:
        doc = {"metrics": {"phonons": {"kappa_103": kappa_103}}}
        mbd_data.round_trip_yaml.dump(doc, file)
    return cast("Model", SimpleNamespace(yaml_path=yaml_path, name="test_model"))


def read_kappa_103(model: Model) -> dict[str, Any]:
    """Read back a model YAML's metrics.phonons.kappa_103."""
    with open(model.yaml_path, encoding="utf-8") as file:
        return mbd_data.round_trip_yaml.load(file)["metrics"]["phonons"]["kappa_103"]


KAPPA_YAML_METRICS = {
    "κ_SRME": 0.25,
    "κ_SRE": 0.125,
    "κ_SRD": -0.05,
    "κ_failure_rate": 0.1,
    "imaginary_mode_rate": 0.05,
    "spectrum_w1": 0.2,
}
AUTHOR_COST = {
    "cost": {"source": "author_reported", "hardware": "A100", "run_time_sec": 28459}
}
OLD_SIDECARS = {
    "force_file": {"name": KAPPA_FORCE, "url": "https://example.com/old-forces.gz"},
    "run_info_file": {"name": KAPPA_RUN_INFO, "url": "https://example.com/old.json"},
}
NEW_PRED_URL = "https://figshare.com/files/new-prediction"
# flat artifact run metadata and the YAML cost mapping it yields
RUN_METADATA = {"cost_source": "runner", "hardware": "H200", "run_time_sec": 12.5}
RUN_COST = {"cost": {"source": "runner", "hardware": "H200", "run_time_sec": 12.5}}


@pytest.mark.parametrize(
    ("initial_metrics", "update_kwargs", "expected"),
    [
        pytest.param(
            {
                "pred_file": {
                    "name": KAPPA_PRED_EXISTING,
                    "url": "https://figshare.com/files/existing",
                },
                "analysis_file": {"name": _GEO_OPT_FILE},
                "force_file": OLD_SIDECARS["force_file"],
            },
            {
                "metrics": KAPPA_METRICS | {"srme": 0.1234, "spectrum_w1": None},
                "force_file_path": KAPPA_FORCE_NEW,
            },
            {
                **KAPPA_YAML_METRICS,
                "κ_SRME": 0.1234,
                "spectrum_w1": None,
                "pred_file": make_file_ref(
                    KAPPA_PRED_EXISTING, url="https://figshare.com/files/existing"
                ),
                "analysis_file": make_file_ref(_GEO_OPT_FILE),
                # a changed sidecar path invalidates its old URL
                "force_file": make_file_ref(KAPPA_FORCE_NEW),
            },
            id="recompute-keeps-artifacts-and-clears-changed-sidecar-url",
        ),
        pytest.param(
            {
                "κ_SRME": 1.0,
                "pred_file": {"name": KAPPA_PRED, "url": "https://example.com/old.gz"},
                **OLD_SIDECARS,
            },
            {
                "run_metadata": RUN_METADATA | {"pred_file_url": NEW_PRED_URL},
                "force_file_path": KAPPA_FORCE,
                "run_info_path": KAPPA_RUN_INFO,
                "replace_pred_file": True,
            },
            {
                **KAPPA_YAML_METRICS,
                **RUN_COST,
                # replaced runs take the run's URL and drop old sidecar URLs even
                # when their paths are unchanged
                "pred_file": make_file_ref(KAPPA_PRED, url=NEW_PRED_URL),
                "force_file": make_file_ref(KAPPA_FORCE),
                "run_info_file": make_file_ref(KAPPA_RUN_INFO),
            },
            id="replace-round-trips-provenance",
        ),
        pytest.param(
            {"κ_SRME": 1.0, **OLD_SIDECARS, **RUN_COST},
            {"replace_pred_file": True},
            {**KAPPA_YAML_METRICS, "pred_file": make_file_ref(KAPPA_PRED)},
            id="replace-clears-stale-sidecars-and-cost",
        ),
        pytest.param(
            # BAM-style author-reported kappa cost survives metric recomputes and
            # new runs whose run_metadata records no cost
            {"κ_SRME": 1.0, "pred_file": {"name": KAPPA_PRED}, **AUTHOR_COST},
            {"run_metadata": {}},
            {**KAPPA_YAML_METRICS, **AUTHOR_COST, "pred_file": {"name": KAPPA_PRED}},
            id="author-reported-cost-kept",
        ),
    ],
)
def test_write_metrics_to_yaml(
    tmp_path: Path,
    initial_metrics: dict[str, object],
    update_kwargs: dict[str, Any],
    expected: dict[str, object],
) -> None:
    """Kappa YAML writes preserve, replace and invalidate artifact metadata."""
    model = kappa_yaml_model(tmp_path, initial_metrics)
    kwargs = dict(update_kwargs)
    metrics = kwargs.pop("metrics", KAPPA_METRICS)
    phonon_metrics.write_metrics_to_yaml(model, metrics, KAPPA_PRED, **kwargs)
    assert read_kappa_103(model) == expected


def test_write_harmonic_to_yaml(tmp_path: Path) -> None:
    """Harmonic-only reruns land in kappa_103.harmonic with their own cost/refs."""
    phonon_path, run_info_path, force_path = (
        f"{_KAPPA_MODEL_DIR}/harmonic/2026-09-05-phonons-kappa-103-{name}"
        for name in ("phonons.json.gz", "run-info.json", "forces.json.gz")
    )
    conductivity_run = {
        "pred_file": {"name": KAPPA_PRED, "url": "https://example.com/pred.gz"},
        **AUTHOR_COST,  # the conductivity run's own cost
        **KAPPA_YAML_METRICS,
    }
    initial = conductivity_run | {
        "phonon_file": {"name": phonon_path},
        "force_file": OLD_SIDECARS["force_file"],
    }
    model = kappa_yaml_model(tmp_path, initial)

    write = partial(
        phonon_metrics.write_harmonic_to_yaml,
        model,
        phonon_file_path=phonon_path,
        run_info_path=run_info_path,
    )
    refs = {
        "phonon_file": {"name": phonon_path},
        "run_info_file": {"name": run_info_path},
    }
    harmonic = write(
        run_metadata=RUN_METADATA | {"n_structures": 103}, force_file_path=force_path
    )
    assert harmonic == refs | {"force_file": {"name": force_path}, **RUN_COST}
    assert read_kappa_103(model) == conductivity_run | {"harmonic": harmonic}

    # a same-name rerun is new content: old url/md5 and cost (tag) never carry over
    uploaded = {"name": phonon_path, "url": "https://example.com/old", "md5": "0" * 32}
    mbd_data.update_yaml_file(
        model.yaml_path,
        "metrics.phonons.kappa_103.harmonic",
        {"phonon_file": uploaded, "cost": {"source": "salvaged"}},
    )
    rerun = RUN_METADATA | {"run_time_sec": 99.0}
    assert write(run_metadata=rerun) == refs | {
        "cost": {**RUN_COST["cost"], "run_time_sec": 99.0}
    }
    with pytest.raises(ValueError, match="records no cost"):
        write(run_metadata={"n_structures": 103})
