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
from matbench_discovery.enums import MbdKey, Model
from matbench_discovery.metrics import phonons as phonon_metrics


@pytest.fixture
def df_pred() -> pd.DataFrame:
    """Mock DataFrame with ML predictions."""
    data = {
        MbdKey.kappa_tot_rta: [np.diag([1, 2, 3]), 2 * np.eye(3)],
        MbdKey.kappa_tot_avg: [2, 2],  # average of diagonal elements
        MbdKey.mode_kappa_tot_rta: [np.diag([0.5, 1, 1.5]), np.eye(3)],
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
        MbdKey.mode_kappa_tot_rta: [np.diag([0.5, 1, 1.5]), np.eye(3)],
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
    ],
)
def test_calculate_kappa_avg_parametrized(
    tensor: NDArray[np.float64], expected: NDArray[np.float64], description: str
) -> None:
    """Test calculation of average thermal conductivity with various inputs."""
    avg = phonon_metrics.calculate_kappa_avg(tensor)
    assert_allclose(avg, expected, rtol=1e-6, err_msg=description)


def test_calculate_kappa_avg_scalar_input() -> None:
    """Scalar input (NaN from a failed calculation) returns NaN, not IndexError."""
    result = phonon_metrics.calculate_kappa_avg(np.asarray(None, dtype=float))
    assert np.all(np.isnan(result))


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
    "n_qpoints,n_modes,weights,expected_len",
    [(4, 3, [1.0, 2.0, 1.0, 1.0], 12), (4, 3, None, None), (4, 3, [1.0, 2.0], None)],
)
def test_mode_weights_for_freqs(
    n_qpoints: int, n_modes: int, weights: list[float] | None, expected_len: int | None
) -> None:
    """q-point weights repeat per mode iff their length matches the q-point axis."""
    ph_freqs = np.ones((n_qpoints, n_modes))
    q_weights = None if weights is None else np.array(weights)
    result = phonon_metrics.mode_weights_for_freqs(ph_freqs, q_weights)
    if expected_len is None or weights is None:
        assert result is None
    else:
        assert result is not None
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

    kappa_srmes = phonon_metrics.calc_kappa_srme(ml_data, dft_data)
    assert kappa_srmes[0] == 0.0  # Should be 0 for identical data

    # Test with different values
    ml_data[MbdKey.kappa_tot_avg] = np.array([3.0])
    ml_data[MbdKey.mode_kappa_tot_rta] = [1.5 * np.eye(3)]
    kappa_srmes = phonon_metrics.calc_kappa_srme(ml_data, dft_data)
    assert kappa_srmes[0] > 0  # Should be non-zero for different values


def test_calc_kappa_metrics_with_different_values(
    df_pred: pd.DataFrame, df_true: pd.DataFrame
) -> None:
    """Test calculation of aggregate metrics with different ML and DFT values."""
    # Modify ML values to be different from DFT
    df_pred_copy = df_pred.copy()
    df_pred_copy[MbdKey.kappa_tot_avg] = [4, 4]  # Double the original values
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
    assert list(phonon_metrics.calc_kappa_srme(ml_data, dft_data)) == [0, 0, 0]
    # Make DFT values half of ML predictions without mutating shared arrays.
    dft_data[MbdKey.kappa_tot_avg] = ml_data[MbdKey.kappa_tot_avg] / 2
    dft_data[MbdKey.kappa_tot_rta] = ml_data[MbdKey.kappa_tot_rta] / 2
    dft_data[MbdKey.mode_kappa_tot_rta] = ml_data[MbdKey.mode_kappa_tot_rta] / 2

    kappa_srmes = phonon_metrics.calc_kappa_srme(ml_data, dft_data)
    assert len(kappa_srmes) == len(ml_data[MbdKey.kappa_tot_avg])
    assert list(kappa_srmes) == pytest.approx([2 / 3, 2 / 3, 2 / 3])

    # The dataframe wrapper keeps one scalar SRME per material, using the first
    # temperature for consistency with SRD/SRE and YAML metrics.
    ml_data[Key.has_imag_ph_modes] = False
    ml_data[Key.final_spg_num] = 1
    ml_data[Key.init_spg_num] = 1

    kappa_srmes = phonon_metrics.calc_kappa_srme_dataframes(
        pd.DataFrame([ml_data]), pd.DataFrame([dft_data])
    )
    assert kappa_srmes == pytest.approx([2 / 3])


def test_calc_kappa_metrics_from_dfs_symmetry(df_minimal: pd.DataFrame) -> None:
    """Test handling of symmetry-related cases in benchmark descriptors."""
    df_pred = pd.concat([df_minimal] * 3, ignore_index=True).copy()
    df_pred[Key.has_imag_ph_modes] = [False, True, False]
    df_pred[Key.final_spg_num] = [1, 1, 2]
    df_pred[Key.init_spg_num] = [1, 1, 1]

    df_true = pd.concat([df_minimal] * 3, ignore_index=True).copy()
    df_true[Key.spg_num] = [1, 1, 1]

    result = phonon_metrics.calc_kappa_metrics_from_dfs(df_pred, df_true)
    assert result[Key.srme].iloc[0] != 2  # Normal case
    assert result[Key.srme].iloc[1] == 2  # Imaginary frequencies
    assert result[Key.srme].iloc[2] == 2  # Broken symmetry
    ignored = phonon_metrics.calc_kappa_metrics_from_dfs(
        df_pred, df_true.copy(), ignore_imaginary_freqs=True
    )
    assert ignored[Key.srme].iloc[1] < 2
    assert bool(df_pred[Key.has_imag_ph_modes].iloc[1]) is True


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
    metrics: dict[str, float] | None = None,
    pred_file_path: str = "models/test/new.json.gz",
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
        {"srme": 0.25, "sre": 0.125} if metrics is None else metrics,
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
            "metrics:\n"
            "  phonons:\n"
            "    kappa_103:\n"
            "      pred_file: models/test/existing.json.gz\n"
            "      pred_file_url: https://figshare.com/files/existing\n"
            "      analysis_file: models/test/analysis.csv.gz\n"
        ),
        metrics={"srme": 0.1234, "sre": 0.5678},
        pred_file_path="models/test/kappa-103.json.gz",
    )
    assert kappa_metrics == {
        "pred_file": "models/test/existing.json.gz",
        "pred_file_url": "https://figshare.com/files/existing",
        "analysis_file": "models/test/analysis.csv.gz",
        "κ_SRME": 0.1234,
        "κ_SRE": 0.5678,
    }


@pytest.mark.parametrize(
    "initial_yaml",
    ["metrics: {}\n", "metrics:\n  phonons: not available\n"],
)
def test_kappa_metric_yaml_replaces_missing_or_unavailable_phonons(
    tmp_path: Path, initial_yaml: str
) -> None:
    """Completed runs replace absent or unavailable phonon metadata."""
    kappa_metrics = update_temp_kappa_yaml(
        tmp_path, initial_yaml, replace_pred_file=True
    )
    assert kappa_metrics["κ_SRME"] == 0.25
    assert kappa_metrics["pred_file"] == "models/test/new.json.gz"


def test_kappa_metric_yaml_round_trip_updates_provenance(tmp_path: Path) -> None:
    """Complete-run metadata and sidecars round-trip through model YAML."""
    kappa_metrics = update_temp_kappa_yaml(
        tmp_path,
        (
            "metrics:\n"
            "  phonons:\n"
            "    kappa_103:\n"
            "      κ_SRME: 1.0\n"
            "      pred_file: models/test/new.json.gz\n"
            "      pred_file_url: https://example.com/old.json.gz\n"
            "      force_file_url: https://example.com/old-forces.json.gz\n"
            "      run_info_file_url: https://example.com/old-run-info.json\n"
        ),
        run_metadata={
            "hardware": "NVIDIA H200",
            "run_time_sec": 12.5,
            "max_rss_gb": 3.5,
            "max_gpu_mem_gb": 4.5,
        },
        force_file_path="models/test/forces.json.gz",
        run_info_path="models/test/run-info.json",
        replace_pred_file=True,
    )
    assert kappa_metrics["κ_SRME"] == 0.25
    assert kappa_metrics["κ_SRE"] == 0.125
    assert kappa_metrics["pred_file"] == "models/test/new.json.gz"
    assert kappa_metrics["pred_file_url"] is None
    assert kappa_metrics["force_file"] == "models/test/forces.json.gz"
    assert kappa_metrics["force_file_url"] is None
    assert kappa_metrics["run_info_file"] == "models/test/run-info.json"
    assert kappa_metrics["run_info_file_url"] is None
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
            "metrics:\n"
            "  phonons:\n"
            "    kappa_103:\n"
            "      κ_SRME: 1.0\n"
            "      pred_file: models/test/pred.json.gz\n"
            "      pred_file_url: https://example.com/pred.json.gz\n"
            "      force_file: models/test/old-forces.json.gz\n"
            "      force_file_url: https://example.com/old-forces.json.gz\n"
        ),
        pred_file_path="models/test/pred.json.gz",
        force_file_path="models/test/new-forces.json.gz",
    )
    assert kappa_metrics["pred_file_url"] == "https://example.com/pred.json.gz"
    assert kappa_metrics["force_file"] == "models/test/new-forces.json.gz"
    assert kappa_metrics["force_file_url"] is None


def test_kappa_metric_yaml_replacement_clears_stale_sidecars(
    tmp_path: Path,
) -> None:
    """Replacing predictions removes provenance absent from the new run."""
    kappa_metrics = update_temp_kappa_yaml(
        tmp_path,
        (
            "metrics:\n"
            "  phonons:\n"
            "    kappa_103:\n"
            "      κ_SRME: 1.0\n"
            "      force_file: models/test/old-forces.json.gz\n"
            "      force_file_url: https://example.com/old-forces.json.gz\n"
            "      run_info_file: models/test/old-run-info.json\n"
            "      run_info_file_url: https://example.com/old-run-info.json\n"
            "      max_gpu_mem_gb: 9.0\n"
        ),
        replace_pred_file=True,
    )
    for stale_key in (
        "force_file",
        "force_file_url",
        "run_info_file",
        "run_info_file_url",
        "max_gpu_mem_gb",
    ):
        assert stale_key not in kappa_metrics
