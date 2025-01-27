"""Tests for thermal conductivity metrics."""

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose
from numpy.typing import NDArray
from pymatviz.enums import Key

from matbench_discovery.enums import MbdKey
from matbench_discovery.metrics import phonons as phonon_metrics


@pytest.fixture
def df_pred() -> pd.DataFrame:
    """Mock DataFrame with ML predictions."""
    data = {
        MbdKey.kappa_tot_rta: [np.diag([1, 2, 3]), 2 * np.eye(3)],
        MbdKey.kappa_tot_avg: [2, 2],  # average of diagonal elements
        MbdKey.mode_kappa_tot: [np.diag([0.5, 1, 1.5]), np.eye(3)],
        Key.mode_weights: [np.array([1]), np.array([1])],
        Key.has_imag_ph_modes: [False, False],
        Key.final_spg_num: [1, 1],
        Key.init_spg_num: [1, 1],
    }
    return pd.DataFrame(data)


@pytest.fixture
def df_true() -> pd.DataFrame:
    """Mock DataFrame with DFT reference values."""
    data = {
        MbdKey.kappa_tot_rta: [np.diag([1, 2, 3]), np.diag([2, 2, 2])],
        MbdKey.kappa_tot_avg: [2, 2],
        MbdKey.mode_kappa_tot: [np.diag([0.5, 1, 1.5]), np.eye(3)],
        Key.mode_weights: [np.array([1]), np.array([1])],
    }
    return pd.DataFrame(data)


@pytest.fixture
def df_minimal() -> pd.DataFrame:
    """Minimal DataFrame with required columns."""
    return pd.DataFrame(
        {
            MbdKey.kappa_tot_rta: [np.eye(3)],
            MbdKey.kappa_tot_avg: [1.0],
            Key.mode_weights: [np.array([1])],
            MbdKey.mode_kappa_tot: [np.eye(3)],
        }
    )


@pytest.fixture
def series_single_temp() -> pd.Series:
    """Mock Series with single temperature data."""
    return pd.Series(
        {
            MbdKey.kappa_tot_avg: np.array([2.0]),
            MbdKey.kappa_tot_rta: 2 * np.eye(3),
            MbdKey.mode_kappa_tot: np.eye(3),
            Key.mode_weights: np.array([1.0]),
        }
    )


@pytest.fixture
def series_multi_temp() -> pd.Series:
    """Mock Series with multi-temperature data."""
    temps = [100, 300, 500]
    return pd.Series(
        {
            MbdKey.kappa_tot_avg: np.array([2.0, 1.5, 1.0]),
            MbdKey.kappa_tot_rta: np.array([2 * np.eye(3), 1.5 * np.eye(3), np.eye(3)]),
            MbdKey.mode_kappa_tot: np.array(
                [np.eye(3), 0.75 * np.eye(3), 0.5 * np.eye(3)]
            ),
            Key.mode_weights: np.ones(len(temps)),
        }
    )


@pytest.mark.parametrize(
    "tensor,expected,description",
    [
        (np.diag([1, 2, 3]), [1 / 3, 2 / 3, 1], "diagonal tensor"),
        (
            [[1, 0.1, 0], [0.1, 2, 0], [0, 0, 3]],
            [0.366667, 0.7, 1],
            "non-diagonal tensor",
        ),
        (np.zeros((3, 3)), np.zeros(3), "zero tensor"),
        (-np.diag([1, 2, 3]), [-1 / 3, -2 / 3, -1], "negative tensor"),
    ],
)
def test_calculate_kappa_avg_parametrized(
    tensor: NDArray[np.float64], expected: NDArray[np.float64], description: str
) -> None:
    """Test calculation of average thermal conductivity with various inputs."""
    avg = phonon_metrics.calculate_kappa_avg(tensor)
    assert_allclose(avg, expected, rtol=1e-6, err_msg=description)


def test_calculate_kappa_avg_edge_cases() -> None:
    """Test average thermal conductivity calculation with edge cases."""
    # Test with zero tensor
    zero_tensor = np.zeros((3, 3))
    result = phonon_metrics.calculate_kappa_avg(zero_tensor)
    assert np.all(result == 0)

    # Test with negative values
    neg_tensor = -np.eye(3)
    result = phonon_metrics.calculate_kappa_avg(neg_tensor)
    assert np.all(result < 0)

    # Test with multiple temperatures
    multi_temp = np.array([np.eye(3), 2 * np.eye(3)])
    result = phonon_metrics.calculate_kappa_avg(multi_temp)
    assert len(result) == 2
    assert np.allclose(result[1], 2 * result[0])

    # Test with NaN values
    tensor_with_nan = np.diag([1.0, 2, 3])  # need dtype=float
    tensor_with_nan[0, 0] = np.nan
    result = phonon_metrics.calculate_kappa_avg(tensor_with_nan)
    assert np.all(np.isnan(result))


@pytest.mark.parametrize(
    "ml_values,dft_values,expected_srd,expected_sre",
    [
        (
            {
                MbdKey.kappa_tot_avg: [1, 2],
                MbdKey.mode_kappa_tot: [np.ones((1, 3, 3))] * 2,
            },
            {
                MbdKey.kappa_tot_avg: [1, 2],
                MbdKey.mode_kappa_tot: [np.ones((1, 3, 3))] * 2,
            },
            [0, 0],
            [0, 0],
        ),
        (
            {
                MbdKey.kappa_tot_avg: [2, 4],
                MbdKey.mode_kappa_tot: [np.ones((1, 3, 3)) * 2] * 2,
            },
            {
                MbdKey.kappa_tot_avg: [1, 2],
                MbdKey.mode_kappa_tot: [np.ones((1, 3, 3))] * 2,
            },
            [2 / 3, 2 / 3],  # (2-1)/1.5, (4-2)/3
            [2 / 3, 2 / 3],
        ),
    ],
)
def test_calc_kappa_metrics_from_dfs_parametrized(
    ml_values: dict[str, list[float]],
    dft_values: dict[str, list[float]],
    expected_srd: list[float],
    expected_sre: list[float],
) -> None:
    """Test processing of benchmark descriptors with various inputs."""
    ml_df = pd.DataFrame(
        {
            **ml_values,
            MbdKey.kappa_tot_rta: [np.ones((3, 3))] * 2,
            Key.mode_weights: [np.array([1])] * 2,
            Key.has_imag_ph_modes: [False] * 2,
            Key.final_spg_num: [1] * 2,
            Key.init_spg_num: [1] * 2,
        }
    )
    dft_df = pd.DataFrame(
        {
            **dft_values,
            MbdKey.kappa_tot_rta: [np.ones((3, 3))] * 2,
            Key.mode_weights: [np.array([1])] * 2,
        }
    )

    result = phonon_metrics.calc_kappa_metrics_from_dfs(ml_df, dft_df)
    assert_allclose(result[Key.srd], expected_srd)
    assert_allclose(result[Key.sre], expected_sre)


@pytest.mark.parametrize(
    "ml_data,expected_srme",
    [
        ({Key.has_imag_ph_modes: True}, 0),  # SRME is 0 for invalid cases
        ({"relaxed_space_group_number": 2}, 0),
        ({MbdKey.kappa_tot_avg: np.array([np.nan])}, [2]),
        (
            {
                MbdKey.kappa_tot_avg: np.array([0]),
                MbdKey.mode_kappa_tot: np.zeros((1, 3, 3)),
            },
            6,  # SRME for zero conductivity case
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
            MbdKey.mode_kappa_tot: np.ones((1, 3, 3)),
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
            MbdKey.mode_kappa_tot: np.ones((1, 3, 3)),
            Key.mode_weights: np.array([1]),
        }
    )

    # Update ml_data with error condition
    for key, val in ml_data.items():
        pred_data[key] = val

    kappa_srmes = phonon_metrics.calc_kappa_srme(pred_data, true_data)
    assert kappa_srmes == pytest.approx(expected_srme)


@pytest.mark.parametrize(
    "temperatures,expected_length",
    [
        ([100.0], 1),
        ([100.0, 200.0], 2),
        ([100.0, 200.0, 300.0], 3),
    ],
)
def test_calc_kappa_srme_temperature_handling_parametrized(
    temperatures: list[float], expected_length: int
) -> None:
    """Test SRME calculation with different numbers of temperatures."""
    ml_data = pd.Series(
        {
            MbdKey.kappa_tot_avg: np.array(temperatures),
            MbdKey.kappa_tot_rta: np.stack([np.eye(3)] * len(temperatures)),
            MbdKey.mode_kappa_tot: np.stack(
                [np.eye(3).reshape(1, 3, 3)] * len(temperatures)
            ),
            Key.mode_weights: np.ones(len(temperatures)),
        }
    )
    dft_data = ml_data.copy()

    kappa_srmes = phonon_metrics.calc_kappa_srme(ml_data, dft_data)
    assert len(kappa_srmes) == expected_length
    assert list(kappa_srmes) == [0] * expected_length


def test_calc_kappa_srme_edge_cases(
    df_pred: pd.DataFrame, df_true: pd.DataFrame
) -> None:
    """Test SRME calculation with various edge cases."""
    # Test imaginary frequencies
    df_imag = df_pred.copy()
    df_imag.loc[0, Key.has_imag_ph_modes] = True
    kappa_srmes = phonon_metrics.calc_kappa_srme_dataframes(df_imag, df_true)
    assert kappa_srmes[0] == 2.0  # Should return 2.0 for imaginary frequencies

    # Test broken symmetry
    df_broken_sym = df_pred.copy()
    df_broken_sym.loc[0, Key.final_spg_num] = 2  # Different from initial
    kappa_srmes = phonon_metrics.calc_kappa_srme_dataframes(df_broken_sym, df_true)
    assert kappa_srmes[0] == 2.0  # Should return 2.0 for broken symmetry

    # Test missing data
    df_missing = df_pred.copy()
    df_missing.loc[0, MbdKey.kappa_tot_avg] = np.nan
    kappa_srmes = phonon_metrics.calc_kappa_srme_dataframes(df_missing, df_true)
    assert kappa_srmes[0] == 2.0  # Should return 2.0 for missing data


def test_calc_kappa_srme_single_material() -> None:
    """Test SRME calculation for a single material."""
    # Create mock data for a single material
    ml_data = pd.Series(
        {
            MbdKey.kappa_tot_avg: np.array([2.0]),
            MbdKey.kappa_tot_rta: 2 * np.eye(3),
            MbdKey.mode_kappa_tot: np.eye(3),
            Key.mode_weights: np.array([1.0]),
        }
    )
    dft_data = ml_data.copy()

    kappa_srmes = phonon_metrics.calc_kappa_srme(ml_data, dft_data)
    assert kappa_srmes[0] == 0.0  # Should be 0 for identical data

    # Test with different values
    ml_data[MbdKey.kappa_tot_avg] = np.array([3.0])
    ml_data[MbdKey.mode_kappa_tot] = [1.5 * np.eye(3)]
    kappa_srmes = phonon_metrics.calc_kappa_srme(ml_data, dft_data)
    assert kappa_srmes[0] > 0  # Should be non-zero for different values


def test_calc_kappa_srme_temperature_handling() -> None:
    """Test SRME calculation with multiple temperatures."""
    # Create mock data with multiple temperatures
    ml_data = pd.Series(
        {
            MbdKey.kappa_tot_avg: np.array([2, 1.5]),  # Two temperatures
            MbdKey.kappa_tot_rta: np.array([2 * np.eye(3), 1.5 * np.eye(3)]),
            MbdKey.mode_kappa_tot: np.array([np.eye(3), 0.75 * np.eye(3)]),
            Key.mode_weights: np.array([1, 1]),
        }
    )
    dft_data = pd.Series(
        {
            MbdKey.kappa_tot_avg: np.array([2, 1.5]),
            MbdKey.kappa_tot_rta: np.array([2 * np.eye(3), 1.5 * np.eye(3)]),
            MbdKey.mode_kappa_tot: np.array([np.eye(3), 0.75 * np.eye(3)]),
            Key.mode_weights: np.array([1, 1]),
        }
    )

    kappa_srmes = phonon_metrics.calc_kappa_srme(ml_data, dft_data)
    assert tuple(kappa_srmes) == (0, 0)  # Should be 0 for identical data


def test_calc_kappa_metrics_with_different_values(
    df_pred: pd.DataFrame, df_true: pd.DataFrame
) -> None:
    """Test calculation of aggregate metrics with different ML and DFT values."""
    # Modify ML values to be different from DFT
    df_pred_copy = df_pred.copy()
    df_pred_copy[MbdKey.kappa_tot_avg] = [4, 4]  # Double the original values
    df_pred_copy[MbdKey.kappa_tot_rta] = [2 * np.diag([1, 2, 3]), 4 * np.eye(3)]
    df_pred_copy[MbdKey.mode_kappa_tot] = [2 * np.diag([1, 2, 3]), 4 * np.eye(3)]

    df_out = phonon_metrics.calc_kappa_metrics_from_dfs(df_pred_copy, df_true)
    assert df_out.shape == (2, 11)
    n_init_cols, n_after_cols = df_pred.shape[1], df_out.shape[1]
    n_new_cols = 4
    assert n_after_cols == n_init_cols + n_new_cols, (
        f"{n_after_cols=} != {n_init_cols=} + {n_new_cols=}"
    )
    pd.testing.assert_index_equal(df_out.index, df_pred.index)
    assert set(df_out) == {
        Key.sre,
        Key.srme,
        Key.final_spg_num,
        Key.has_imag_ph_modes,
        Key.init_spg_num,
        MbdKey.kappa_tot_avg,
        MbdKey.kappa_tot_rta,
        MbdKey.mode_kappa_tot_avg,
        Key.mode_weights,
        Key.srd,
        MbdKey.true_kappa_tot_avg,
    }
    assert df_out[Key.sre].mean() == pytest.approx(2 / 3)
    assert df_out[Key.srme].mean() == pytest.approx(1 / 4)


def test_calc_kappa_metrics_from_dfs_missing_columns(
    df_minimal: pd.DataFrame,
) -> None:
    """Test processing benchmark descriptors with missing columns."""
    df_pred = df_minimal.copy()
    df_pred = df_pred.drop(columns=[MbdKey.kappa_tot_avg])
    df_true = df_minimal.copy()

    df_out = phonon_metrics.calc_kappa_metrics_from_dfs(df_pred, df_true)
    assert MbdKey.kappa_tot_avg in df_out
    assert Key.srd in df_out
    assert Key.sre in df_out
    assert Key.srme in df_out


def test_calc_kappa_srme_temperature_dependence(series_multi_temp: pd.Series) -> None:
    """Test SRME calculation with temperature-dependent conductivities."""
    ml_data = series_multi_temp.copy()
    dft_data = series_multi_temp.copy()
    dft_data[MbdKey.kappa_tot_avg] /= 2  # Make DFT values half of ML predictions
    dft_data[MbdKey.kappa_tot_rta] /= 2
    dft_data[MbdKey.mode_kappa_tot] /= 2

    kappa_srmes = phonon_metrics.calc_kappa_srme(ml_data, dft_data)
    assert len(kappa_srmes) == len(ml_data[Key.mode_weights])
    # TODO Should be non-zero since ML predictions are double DFT
    assert kappa_srmes.tolist() == [0, 0, 0]


def test_calc_kappa_metrics_from_dfs_symmetry(df_minimal: pd.DataFrame) -> None:
    """Test handling of symmetry-related cases in benchmark descriptors."""
    df_pred = pd.concat([df_minimal] * 3, ignore_index=True)
    df_pred[Key.has_imag_ph_modes] = [False, True, False]
    df_pred[Key.final_spg_num] = [1, 1, 2]
    df_pred[Key.init_spg_num] = [1, 1, 1]

    df_true = pd.concat([df_minimal] * 3, ignore_index=True)
    df_true[Key.spg_num] = [1, 1, 1]

    result = phonon_metrics.calc_kappa_metrics_from_dfs(df_pred, df_true)
    assert result[Key.srme].iloc[0] != 2  # Normal case
    assert result[Key.srme].iloc[1] == 2  # Imaginary frequencies
    assert result[Key.srme].iloc[2] == 2  # Broken symmetry


def test_calc_kappa_srme_dataframes_error_handling(df_minimal: pd.DataFrame) -> None:
    """Test error handling in SRME calculation for dataframes."""
    df_pred = pd.concat([df_minimal] * 2, ignore_index=True)
    df_pred.loc[0, MbdKey.kappa_tot_avg] = np.nan
    df_pred[Key.has_imag_ph_modes] = [True, False]
    df_pred[Key.final_spg_num] = [2, 1]
    df_pred[Key.init_spg_num] = [1, 1]

    df_true = pd.concat([df_minimal] * 2, ignore_index=True)
    df_true[Key.spg_num] = [1, 1]

    result = phonon_metrics.calc_kappa_srme_dataframes(df_pred, df_true)
    assert len(result) == 2
    assert result[0] == 2  # First entry should be 2 due to imaginary frequencies
    assert 0 <= result[1] <= 2  # Second entry should be valid SRME value
