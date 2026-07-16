"""Tests for thermal conductivity metrics."""

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
from matbench_discovery.data import file_ref_url, make_file_ref
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
    values = np.array([1.0, 2.0, 4.0, 8.0])
    weights = np.array([1.0, 3.0, 2.0, 1.0])
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


def test_calc_kappa_metrics_from_dfs() -> None:
    """Aggregate metrics compare predicted tensors with reference scalars."""
    ml_df = pd.DataFrame(
        {
            MbdKey.mode_kappa_tot_rta: [np.ones((1, 3, 3))] * 2,
            MbdKey.kappa_tot_rta: [np.ones((3, 3))] * 2,
            Key.mode_weights: [np.array([1])] * 2,
            Key.has_imag_ph_modes: [False] * 2,
            Key.final_spg_num: [1] * 2,
            Key.init_spg_num: [1] * 2,
        }
    )
    dft_df = pd.DataFrame(
        {
            MbdKey.kappa_tot_avg: [1, 2],
            MbdKey.mode_kappa_tot_rta: [np.ones((1, 3, 3))] * 2,
            MbdKey.kappa_tot_rta: [np.ones((3, 3))] * 2,
            Key.mode_weights: [np.array([1])] * 2,
        }
    )

    result = phonon_metrics.calc_kappa_metrics_from_dfs(ml_df, dft_df)
    assert_allclose(result[Key.srd], [0, -2 / 3])
    assert_allclose(result[Key.sre], [0, 2 / 3])


@pytest.mark.parametrize(
    "ml_data,expected_srme",
    [
        ({Key.has_imag_ph_modes: True}, 0),  # SRME is 0 for invalid cases
        ({"relaxed_space_group_number": 2}, 0),
        ({MbdKey.kappa_tot_avg: np.array([np.nan])}, [2]),
        (
            {
                MbdKey.kappa_tot_avg: np.array([0]),
                MbdKey.mode_kappa_tot_rta: np.zeros((1, 3, 3)),
            },
            [2],
        ),
    ],
)
def test_calc_kappa_srme_error_cases(
    ml_data: dict[str, list[float]], expected_srme: float
) -> None:
    """Test SRME calculation with various error conditions."""
    pred_data = pd.Series(
        {
            MbdKey.kappa_tot_avg: np.array([1]),
            MbdKey.kappa_tot_rta: np.ones((3, 3)),
            MbdKey.mode_kappa_tot_rta: np.ones((1, 3, 3)),
            Key.mode_weights: np.array([1]),
            Key.has_imag_ph_modes: False,
            Key.final_spg_num: 1,
            Key.init_spg_num: 1,
        }
    )
    true_data = pd.Series(
        {
            MbdKey.kappa_tot_avg: np.array([1]),
            MbdKey.kappa_tot_rta: np.ones((3, 3)),
            MbdKey.mode_kappa_tot_rta: np.ones((1, 3, 3)),
            Key.mode_weights: np.array([1]),
        }
    )

    # Update ml_data with error condition
    for key, val in ml_data.items():
        pred_data[key] = val

    kappa_srmes = phonon_metrics.calc_kappa_srme(pred_data, true_data)
    assert kappa_srmes == pytest.approx(expected_srme)


@pytest.mark.parametrize(
    ("column", "value"),
    [
        (Key.has_imag_ph_modes, True),
        (Key.final_spg_num, 2),
        (MbdKey.kappa_tot_avg, np.nan),
        (MbdKey.kappa_tot_avg, np.inf),
        (MbdKey.kappa_tot_avg, -1),
    ],
)
def test_calc_kappa_srme_invalid_predictions_score_two(
    df_pred: pd.DataFrame,
    df_true: pd.DataFrame,
    column: str,
    value: object,
) -> None:
    """Imaginary modes, broken symmetry, and missing data receive maximum SRME."""
    modified = df_pred.copy()
    modified[column] = modified[column].astype(object)
    modified.loc[0, column] = value
    assert phonon_metrics.calc_kappa_srme_dataframes(modified, df_true)[0] == 2


def test_calc_kappa_srme_single_material() -> None:
    """Test SRME calculation for a single material."""
    # Create mock data for a single material
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

    # Test with different values
    ml_data[MbdKey.kappa_tot_avg] = np.array([3.0])
    ml_data[MbdKey.mode_kappa_tot_rta] = [1.5 * np.eye(3)]
    assert phonon_metrics.calc_kappa_srme(ml_data, dft_data)[0] > 0


def test_calc_kappa_metrics_with_different_values(
    df_pred: pd.DataFrame, df_true: pd.DataFrame
) -> None:
    """Test calculation of aggregate metrics with different ML and DFT values."""
    # Modify ML values to be different from DFT
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
        MbdKey.true_kappa_tot_avg,
    }
    assert df_out[Key.sre].mean() == pytest.approx(2 / 3)
    assert df_out[Key.srme].mean() == pytest.approx(1)


def test_calc_kappa_metrics_from_dfs_missing_columns(
    df_minimal: pd.DataFrame,
) -> None:
    """Test processing benchmark descriptors with missing columns."""
    df_pred = df_minimal.drop(columns=[MbdKey.kappa_tot_avg])
    df_true = df_minimal.copy()

    df_out = phonon_metrics.calc_kappa_metrics_from_dfs(df_pred, df_true)
    assert {MbdKey.kappa_tot_avg, Key.srd, Key.sre, Key.srme} <= set(df_out)


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
    """Test handling of symmetry-related cases in benchmark descriptors."""
    df_pred = pd.concat([df_minimal] * 3, ignore_index=True).copy()
    df_pred[Key.has_imag_ph_modes] = [False, True, False]
    df_pred[Key.final_spg_num] = [1, 1, 2]
    df_pred[Key.init_spg_num] = [1, 1, 1]

    df_true = pd.concat([df_minimal] * 3, ignore_index=True).copy()
    df_true[Key.spg_num] = [1, 1, 1]

    result = phonon_metrics.calc_kappa_metrics_from_dfs(df_pred, df_true)
    assert result[Key.srme].eq(2).tolist() == [False, True, True]


@pytest.mark.parametrize(
    ("reference_spg_col", "expected_failures"),
    [
        (Key.init_spg_num, [False, True]),
        (Key.spg_num, [False, True]),
        (None, [False, False]),
    ],
)
def test_calc_kappa_srme_dataframes_missing_init_spg_uses_reference(
    df_minimal: pd.DataFrame,
    reference_spg_col: str | None,
    expected_failures: list[bool],
) -> None:
    """Missing prediction symmetry uses canonical, legacy, or absent references."""
    df_pred = pd.concat([df_minimal] * 2, ignore_index=True).copy()
    df_pred[Key.final_spg_num] = [225, 186]
    df_true = pd.concat([df_minimal] * 2, ignore_index=True).copy()
    if reference_spg_col is not None:
        df_true[reference_spg_col] = [225, 225]

    srme_values = phonon_metrics.calc_kappa_srme_dataframes(df_pred, df_true)
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

    srme_values = phonon_metrics.calc_kappa_srme_dataframes(df_pred, df_true)
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


def test_evaluate_kappa_predictions_uses_exact_failure_boundary(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Near-limit SRME is not a failure and an absent imaginary flag counts as zero."""
    from matbench_discovery import phonons

    material_id = "mp-1"
    reference = pd.DataFrame(index=[material_id])
    predictions = pd.DataFrame(index=[material_id])

    def read_reference(_path: str) -> pd.DataFrame:
        """Return the fixed one-material reference frame."""
        return reference

    def calculate_metrics(
        _predictions: pd.DataFrame, _reference: pd.DataFrame
    ) -> pd.DataFrame:
        """Return a valid score immediately below the hard failure boundary."""
        return pd.DataFrame(
            {Key.srme: [1.9999], Key.sre: [1.9999], Key.srd: [1.9999]},
            index=[material_id],
        )

    monkeypatch.setattr(phonons, "read_kappa_json", read_reference)
    monkeypatch.setattr(
        phonon_metrics, "calc_kappa_metrics_from_dfs", calculate_metrics
    )
    metrics = phonon_metrics.evaluate_kappa_predictions(predictions)

    assert metrics["failure_rate"] == 0
    assert metrics["imaginary_mode_rate"] == 0


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


def test_kappa_analysis_preserves_hard_failure_boundary() -> None:
    """Rounded payload SRME equals 2 only for exactly censored failures."""
    from scripts.model_figs.kappa_103_analysis import first_scalar

    assert first_scalar(2.0, exclusive_upper_bound=2.0) == 2
    assert first_scalar(1.99996, exclusive_upper_bound=2.0) == 1.9999


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
    df_pred: pd.DataFrame,
    df_true: pd.DataFrame,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The analysis payload delegates spectrum W1 to the canonical metric helper."""
    from scripts.model_figs.kappa_103_analysis import model_payload

    material_ids = ["mp-1", "mp-2"]
    df_pred = df_pred.set_axis(material_ids)
    df_true = df_true.set_axis(material_ids)
    frequencies = [np.array([[0.0, 1.0, 2.0]])] * len(material_ids)
    df_pred[Key.ph_freqs] = frequencies
    df_true[Key.ph_freqs] = frequencies
    call_count = 0

    def shared_w1(_prediction: pd.Series, _reference: pd.Series) -> float:
        """Return a sentinel W1 value while counting delegated calls."""
        nonlocal call_count
        call_count += 1
        return 0.12345

    monkeypatch.setattr(phonon_metrics, "phonon_spectrum_w1", shared_w1)
    payload = model_payload(Model.mace_mp_0, df_pred, df_true, material_ids)

    assert call_count == len(material_ids)
    assert payload["freq_w1"] == [0.1235, 0.1235]
    assert payload["freq_w1_mean"] == 0.1235


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


def update_temp_kappa_yaml(
    tmp_path: Path,
    initial_yaml: str,
    *,
    metrics: dict[str, float | None] | None = None,
    pred_file_path: str = KAPPA_PRED,
    run_metadata: dict[str, object] | None = None,
    force_file_path: str | None = None,
    run_info_path: str | None = None,
    replace_pred_file: bool = False,
) -> dict[str, Any]:
    """Update a temporary model YAML and return its kappa metrics mapping."""
    yaml_path = f"{tmp_path}/model.yml"
    with open(yaml_path, mode="w", encoding="utf-8") as file:
        file.write(initial_yaml)
    model = cast("Model", SimpleNamespace(yaml_path=yaml_path))
    phonon_metrics.write_metrics_to_yaml(
        model,
        KAPPA_METRICS if metrics is None else metrics,
        pred_file_path,
        run_metadata=run_metadata,
        force_file_path=force_file_path,
        run_info_path=run_info_path,
        replace_pred_file=replace_pred_file,
    )
    with open(yaml_path, encoding="utf-8") as file:
        metadata = mbd_data.round_trip_yaml.load(file)
    return cast("dict[str, Any]", metadata["metrics"]["phonons"]["kappa_103"])


def test_write_metrics_to_yaml_preserves_existing_artifacts(tmp_path: Path) -> None:
    """Metric recomputation preserves established artifact metadata."""
    kappa_metrics = update_temp_kappa_yaml(
        tmp_path,
        (
            f"metrics:\n  phonons:\n    kappa_103:\n"
            f"      pred_file:\n"
            f"        name: {KAPPA_PRED_EXISTING}\n"
            f"        url: https://figshare.com/files/existing\n"
            f"      analysis_file:\n"
            f"        name: {_GEO_OPT_FILE}\n"
        ),
        metrics={
            **KAPPA_METRICS,
            "srme": 0.1234,
            "sre": 0.5678,
            "spectrum_w1": None,
        },
        pred_file_path=KAPPA_PRED,
    )
    assert kappa_metrics == {
        "pred_file": make_file_ref(
            KAPPA_PRED_EXISTING, url="https://figshare.com/files/existing"
        ),
        "analysis_file": make_file_ref(_GEO_OPT_FILE),
        "κ_SRME": 0.1234,
        "κ_SRE": 0.5678,
        "κ_SRD": -0.05,
        "κ_failure_rate": 0.1,
        "imaginary_mode_rate": 0.05,
        "spectrum_w1": None,
    }


@pytest.mark.parametrize(
    "pred_file_url",
    [None, "https://figshare.com/files/new-prediction"],
    ids=["without-url", "with-run-url"],
)
def test_kappa_metric_yaml_creates_prediction_file(
    tmp_path: Path, pred_file_url: str | None
) -> None:
    """Completed runs create phonon metadata with optional authoritative URLs."""
    kappa_metrics = update_temp_kappa_yaml(
        tmp_path,
        "metrics: {}\n",
        run_metadata=(
            {"pred_file_url": pred_file_url} if pred_file_url is not None else None
        ),
        replace_pred_file=True,
    )
    assert kappa_metrics["κ_SRME"] == 0.25
    assert kappa_metrics["pred_file"] == make_file_ref(KAPPA_PRED, url=pred_file_url)


def test_kappa_metric_yaml_round_trip_updates_provenance(tmp_path: Path) -> None:
    """Complete-run metadata and sidecars round-trip through model YAML."""
    kappa_metrics = update_temp_kappa_yaml(
        tmp_path,
        (
            f"metrics:\n  phonons:\n    kappa_103:\n"
            f"      κ_SRME: 1.0\n"
            f"      pred_file:\n"
            f"        name: {KAPPA_PRED}\n"
            f"        url: https://example.com/old.json.gz\n"
            f"      force_file:\n"
            f"        name: {KAPPA_FORCE}\n"
            f"        url: https://example.com/old-forces.json.gz\n"
            f"      run_info_file:\n"
            f"        name: {KAPPA_RUN_INFO}\n"
            f"        url: https://example.com/old-run-info.json\n"
        ),
        run_metadata={
            "hardware": "NVIDIA H200",
            "run_time_sec": 12.5,
            "max_rss_gb": 3.5,
            "max_gpu_mem_gb": 4.5,
        },
        force_file_path=KAPPA_FORCE,
        run_info_path=KAPPA_RUN_INFO,
        replace_pred_file=True,
    )
    assert kappa_metrics["κ_SRME"] == 0.25
    assert kappa_metrics["κ_SRE"] == 0.125
    assert kappa_metrics["pred_file"] == make_file_ref(KAPPA_PRED)
    assert kappa_metrics["force_file"] == make_file_ref(KAPPA_FORCE)
    assert kappa_metrics["run_info_file"] == make_file_ref(KAPPA_RUN_INFO)
    assert kappa_metrics["hardware"] == "NVIDIA H200"
    assert kappa_metrics["run_time_sec"] == 12.5
    assert kappa_metrics["max_rss_gb"] == 3.5
    assert kappa_metrics["max_gpu_mem_gb"] == 4.5


def test_kappa_metric_yaml_clears_url_when_sidecar_path_changes(
    tmp_path: Path,
) -> None:
    """Changing a sidecar path invalidates its existing remote URL."""
    kappa_metrics = update_temp_kappa_yaml(
        tmp_path,
        (
            f"metrics:\n  phonons:\n    kappa_103:\n"
            f"      κ_SRME: 1.0\n"
            f"      pred_file:\n"
            f"        name: {KAPPA_PRED}\n"
            f"        url: https://example.com/pred.json.gz\n"
            f"      force_file:\n"
            f"        name: {KAPPA_FORCE}\n"
            f"        url: https://example.com/old-forces.json.gz\n"
        ),
        pred_file_path=KAPPA_PRED,
        force_file_path=KAPPA_FORCE_NEW,
    )
    assert (
        file_ref_url(kappa_metrics["pred_file"]) == "https://example.com/pred.json.gz"
    )
    assert kappa_metrics["force_file"] == make_file_ref(KAPPA_FORCE_NEW)


def test_kappa_metric_yaml_replacement_clears_stale_sidecars(
    tmp_path: Path,
) -> None:
    """Replacing predictions removes provenance absent from the new run."""
    kappa_metrics = update_temp_kappa_yaml(
        tmp_path,
        (
            f"metrics:\n  phonons:\n    kappa_103:\n"
            f"      κ_SRME: 1.0\n"
            f"      force_file:\n"
            f"        name: {KAPPA_FORCE}\n"
            f"        url: https://example.com/old-forces.json.gz\n"
            f"      run_info_file:\n"
            f"        name: {KAPPA_RUN_INFO}\n"
            f"        url: https://example.com/old-run-info.json\n"
            f"      max_gpu_mem_gb: 9.0\n"
        ),
        replace_pred_file=True,
    )
    for stale_key in ("force_file", "run_info_file", "max_gpu_mem_gb"):
        assert stale_key not in kappa_metrics
