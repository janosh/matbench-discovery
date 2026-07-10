"""Test diatomic curve metrics calculation functions."""

import re
from collections.abc import Callable

import numpy as np
import pytest

from matbench_discovery.metrics import diatomics
from matbench_discovery.metrics.diatomics import DiatomicCurves
from matbench_discovery.metrics.diatomics.energy import (
    _threshold_diff_signs,
    _validate_diatomic_curve,
    calc_energy_diff_flips,
    calc_energy_jump,
    calc_pbe_bond_length_error,
    calc_pbe_energy_mae,
    calc_pbe_vib_freq_error,
    calc_pbe_wall_dist_mae,
    calc_pbe_well_depth_error,
)


@pytest.mark.parametrize(
    "metric_func",
    [
        diatomics.calc_tortuosity,
        diatomics.calc_energy_jump,
        diatomics.calc_energy_diff_flips,
    ],
)
def test_energy_metrics_on_fixture(
    metric_func: Callable[..., float],
    pred_ref_diatomic_curves: tuple[DiatomicCurves, DiatomicCurves],
) -> None:
    """All single-curve energy metrics return non-negative floats on fixture data."""
    _ref_curves, pred_curves = pred_ref_diatomic_curves
    pred_curve = pred_curves.homo_nuclear["H"]
    seps_pred, energies_pred = pred_curve.distances, pred_curve.energies
    result = metric_func(seps_pred, energies_pred)
    assert isinstance(result, float)
    assert result >= 0


@pytest.mark.parametrize(
    "vals, threshold, expected_n_flips, expected_diffs_len",
    [
        # Monotone increasing: no flips
        (np.array([1.0, 2.0, 3.0, 4.0]), 1e-3, 0, 3),
        # Single oscillation: 1 flip (up then down)
        (np.array([1.0, 3.0, 2.0]), 1e-3, 1, 2),
        # Multiple oscillations: 3 flips
        (np.array([1.0, 3.0, 2.0, 4.0, 1.0]), 1e-3, 3, 4),
        # Diffs below threshold get zeroed → no flips
        (np.array([1.0, 1.0004, 1.0002, 1.0006]), 1e-3, 0, 0),
        # Two points: 1 diff, 0 flips
        (np.array([1.0, 2.0]), 1e-3, 0, 1),
    ],
)
def test_threshold_diff_signs(
    vals: np.ndarray, threshold: float, expected_n_flips: int, expected_diffs_len: int
) -> None:
    """Test _threshold_diff_signs with known inputs."""
    diffs, signs, flips = _threshold_diff_signs(vals, threshold)
    assert len(diffs) == expected_diffs_len
    assert len(signs) == expected_diffs_len
    assert int(np.sum(flips)) == expected_n_flips
    # signs should only be -1 or +1 (zeros are filtered)
    if len(signs) > 0:
        assert set(signs.tolist()) <= {-1.0, 1.0}


@pytest.mark.parametrize(
    "energies, expected_flips, expected_jump",
    [
        # Monotone: 0 flips, 0 jump
        (np.array([1.0, 2.0, 3.0, 4.0, 5.0]), 0, 0.0),
        # diffs=[2,-1,2,1], signs=[+,-,+,+], flips at positions 0,1 → 2 flips
        # jump = |diffs[0]|+|diffs[1]| + |diffs[1]|+|diffs[2]| = 2+1+1+2 = 6
        (np.array([1.0, 3.0, 2.0, 4.0, 5.0]), 2, 6.0),
        # diffs=[2,-1,2,-2.5], signs=[+,-,+,-], 3 flips
        # jump = |2|+|1| + |1|+|2| + |2|+|2.5| = 3+3+4.5 = 10.5
        (np.array([0.0, 2.0, 1.0, 3.0, 0.5]), 3, 10.5),
    ],
)
def test_energy_flips_and_jumps_concrete(
    energies: np.ndarray, expected_flips: int, expected_jump: float
) -> None:
    """Test calc_energy_diff_flips and calc_energy_jump with concrete values."""
    seps = np.arange(1, len(energies) + 1, dtype=float)
    assert calc_energy_diff_flips(seps, energies) == expected_flips
    assert calc_energy_jump(seps, energies) == pytest.approx(expected_jump)


def test_validate_normalize_energy() -> None:
    """Test that normalize_energy shifts energies to zero at far field."""
    seps = np.array([1.0, 2.0, 3.0, 4.0])
    energies = np.array([10.0, 5.0, 2.0, 1.0])
    _, normed = _validate_diatomic_curve(seps, energies, normalize_energy=True)
    # After sort ascending, last value (at sep=4) should be 0
    assert normed[-1] == 0.0
    assert normed[0] == pytest.approx(9.0)  # 10 - 1

    # Forces (ndim > 1) should NOT be normalized even with flag set
    forces = np.ones((4, 2, 3))
    _, forces_out = _validate_diatomic_curve(seps, forces, normalize_energy=True)
    np.testing.assert_array_equal(forces_out, forces)


def test_pbe_reference_energy_metrics() -> None:
    """PBE-reference metrics on parabolic wells with analytically-known parameters.

    ref = (r - r_ref_eq)^2 + E_ref_min, pred = k_pred * (r - r_pred_eq)^2 + E_pred_min
    on 0.05 A grids that hit both minima exactly, so bond-length and well-depth errors
    have closed forms. Wall-dist/energy-MAE/vib-freq expectations are regression pins
    (they depend on interpolation grids and physical constants).
    """
    r_ref_eq, e_ref_min, r_ref_max = 1.5, -2, 3.0
    k_pred, r_pred_eq, e_pred_min, r_pred_max = 1.2, 1.6, -1.7, 3.05
    seps_ref = np.linspace(0.5, r_ref_max, 51)
    energies_ref = (seps_ref - r_ref_eq) ** 2 + e_ref_min
    seps_pred = np.linspace(0.55, r_pred_max, 51)
    energies_pred = k_pred * (seps_pred - r_pred_eq) ** 2 + e_pred_min

    assert calc_pbe_wall_dist_mae(
        seps_ref, energies_ref, seps_pred, energies_pred, thresholds_ev=(1,)
    ) == pytest.approx(0.187, abs=0.01)
    assert calc_pbe_energy_mae(
        seps_ref, energies_ref, seps_pred, energies_pred, interpolate=200
    ) == pytest.approx(0.118, abs=0.01)
    assert calc_pbe_bond_length_error(
        seps_ref, energies_ref, seps_pred, energies_pred
    ) == pytest.approx(r_pred_eq - r_ref_eq)
    # well depth D_e = E(r_max) - E_min, so the error is the curvature-scaled
    # difference of squared distances from equilibrium to the last grid point
    expected_well_depth_err = abs(
        k_pred * (r_pred_max - r_pred_eq) ** 2 - (r_ref_max - r_ref_eq) ** 2
    )
    assert calc_pbe_well_depth_error(
        seps_ref, energies_ref, seps_pred, energies_pred
    ) == pytest.approx(expected_well_depth_err, abs=0.01)
    assert calc_pbe_vib_freq_error(
        "H", seps_ref, energies_ref, seps_pred, energies_pred
    ) == pytest.approx(99.1, abs=1)


def test_wall_dist_penalizes_unreached_prediction_threshold() -> None:
    """A too-soft predicted wall receives the full reference-radius error."""
    separations = np.array([0.5, 1.0, 1.5, 2.0])
    reference_energies = np.array([100.0, 25.0, 0.0, 0.0])
    predicted_energies = np.array([10.0, 5.0, 0.0, 0.0])

    assert calc_pbe_wall_dist_mae(
        separations,
        reference_energies,
        separations,
        predicted_energies,
        thresholds_ev=(50,),
    ) == pytest.approx(5 / 6)


def test_pbe_reference_metrics_skip_unbound_refs() -> None:
    """Well-based PBE metrics skip unbound reference curves."""
    seps = np.linspace(1, 5, 20)
    flat_energies = np.zeros_like(seps)
    pred_energies = (seps - 2) ** 2
    for metric_func in (calc_pbe_bond_length_error, calc_pbe_well_depth_error):
        assert np.isnan(metric_func(seps, flat_energies, seps, pred_energies))
    assert np.isnan(
        calc_pbe_vib_freq_error("He", seps, flat_energies, seps, pred_energies)
    )


@pytest.mark.parametrize(
    "energies, expected",
    [
        (np.arange(1, 6), 1.0),
        (np.arange(1, 6) ** 2, 1.0),
        (np.ones(5), np.nan),
        (np.array([1, 1, 1.5, 2, 2]), 1.0),
    ],
    ids=["linear", "quadratic", "constant", "step"],
)
def test_tortuosity_simple_cases(energies: np.ndarray, expected: float) -> None:
    """Test tortuosity with simple, easy-to-reason-about cases."""
    dists = np.arange(1, 6)
    result = diatomics.calc_tortuosity(dists, energies)
    assert result == pytest.approx(expected, nan_ok=True)


@pytest.mark.parametrize(
    "xs, ys, match",
    [
        (np.array([1, np.nan, 3, 4, 5]), np.arange(1, 6), "Input contains NaN"),
        (np.arange(1, 6), np.array([1, np.inf, 3, 4, 5]), "Input contains infinite"),
        (np.array([1]), np.array([1]), "Input must have at least 2 points"),
        (np.array([1, 1, 2, 3, 4]), np.arange(1, 6), "xs contains 1 duplicates"),
        # non-adjacent duplicates (regression: np.diff used to miss these)
        (np.array([1, 3, 2, 1, 4]), np.arange(1, 6), "xs contains 1 duplicates"),
        (np.array([1, 2, 1, 3, 2]), np.arange(1, 6), "xs contains 2 duplicates"),
        (
            np.array([1, 2, 3]),
            np.arange(1, 6),
            re.escape("len(xs_arr)=3 != len(ys_arr)=5"),
        ),
    ],
)
def test_validation_errors(xs: np.ndarray, ys: np.ndarray, match: str) -> None:
    """Test that _validate_diatomic_curve rejects invalid inputs."""
    with pytest.raises(ValueError, match=match):
        calc_energy_jump(xs, ys)


def test_unsorted_input_accepted() -> None:
    """Unsorted x values should work since sorting is handled internally."""
    xs = np.array([1, 3, 2, 5, 4])
    ys = np.array([1, 3, 2, 5, 4])
    assert np.isfinite(calc_pbe_energy_mae(np.arange(1, 6), np.arange(1, 6), xs, ys))
