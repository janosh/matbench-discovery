"""Test diatomic curve metrics calculation functions."""

import re
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from matbench_discovery.enums import MbdKey, Model
from matbench_discovery.metrics import diatomics
from matbench_discovery.metrics.diatomics import DiatomicCurve, DiatomicCurves

np_rng = np.random.default_rng(seed=0)


def test_diatomic_classes() -> None:
    """Test DiatomicCurve and DiatomicCurves initialization and validation."""
    dists = np.linspace(0.5, 5.0, 10).tolist()
    energies = np_rng.random(len(dists)).tolist()
    forces = np_rng.random((len(dists), 2, 3)).tolist()

    # Test DiatomicCurve initialization and array conversion
    curve = DiatomicCurve(distances=dists, energies=energies, forces=forces)
    assert {*map(type, (curve.distances, curve.energies, curve.forces))} == {np.ndarray}
    for orig, processed in [
        (curve.distances, dists),
        (curve.energies, energies),
        (curve.forces, forces),
    ]:
        np.testing.assert_array_equal(orig, processed)

    # Test DiatomicCurves initialization and from_dict
    data = {
        "distances": dists,
        "homo-nuclear": {"H": {"energies": energies, "forces": forces}},
        "hetero-nuclear": {"H-He": {"energies": energies, "forces": forces}},
    }
    curves = DiatomicCurves.from_dict(data)
    assert isinstance(curves.homo_nuclear["H"], DiatomicCurve)
    print(f"{curves.hetero_nuclear=}")
    h_he_curve = curves.hetero_nuclear.get("H-He")
    assert isinstance(h_he_curve, DiatomicCurve)
    np.testing.assert_array_equal(curves.distances, dists)
    np.testing.assert_array_equal(curves.homo_nuclear["H"].energies, energies)
    np.testing.assert_array_equal(h_he_curve.energies, energies)


def base_curve(xs: np.ndarray) -> np.ndarray:
    """Simple Morse potential."""
    return 5 * (1 - np.exp(-2 * (xs - 1.5))) ** 2 - 5


@pytest.mark.parametrize(
    "expected",
    [
        {
            "name": "no_mod",
            MbdKey.norm_auc: 0.0,
            MbdKey.energy_mae: 0.0,
            MbdKey.tortuosity: 1.0,
            MbdKey.energy_jump: 0.083043538,
            MbdKey.smoothness: 578.772692,
        },
        {
            "name": "vertical_shift",
            MbdKey.norm_auc: 0.00489955,
            MbdKey.norm_auc: 0.00489955,
            MbdKey.energy_mae: 1.0,
            MbdKey.tortuosity: 1.0,
            MbdKey.energy_jump: 0.083043538,
            MbdKey.smoothness: 578.772692,
        },
        {
            "name": "morse",
            MbdKey.norm_auc: 0.411922297,
            MbdKey.energy_mae: 84.797947,
            MbdKey.tortuosity: 1.0,
            MbdKey.energy_jump: 0.083043538,
            MbdKey.smoothness: 4561.131837,
        },
    ],
)
def test_curve_shifts(expected: dict[str, float | str]) -> None:
    """Test metrics for various curve modifications."""
    name = str(expected.pop("name"))
    curve_func = {
        "no_mod": lambda xs: base_curve(xs),
        "vertical_shift": lambda xs: base_curve(xs) + 1.0,
        "morse": lambda xs: 5 * (1 - np.exp(-2 * (xs - 2.0))) ** 2 - 5,
    }[name]
    dists = np.linspace(0.5, 5.0, 100)
    e_ref = base_curve(dists)
    e_pred = curve_func(dists)

    # Create forces array with shape (n_distances, n_atoms, 3)
    dummy_forces = np.zeros((len(dists), 2, 3))  # dummy forces, unused
    dummy_forces[:, 0, 0] = -np.gradient(e_ref, dists)  # force on first atom
    dummy_forces[:, 1, 0] = np.gradient(e_ref, dists)  # force on second atom

    h_ref = DiatomicCurve(distances=dists, energies=e_ref, forces=dummy_forces)
    ref_curves = DiatomicCurves(distances=dists, homo_nuclear={"H": h_ref})

    h_pred = DiatomicCurve(distances=dists, energies=e_pred, forces=dummy_forces)
    pred_curves = DiatomicCurves(distances=dists, homo_nuclear={"H": h_pred})

    metrics_out = diatomics.calc_diatomic_metrics(ref_curves, pred_curves)

    # Check metrics
    for metric_key, expect in expected.items():
        actual = metrics_out["H"][metric_key]
        assert actual == pytest.approx(expect), (
            f"{name=}, {metric_key=} expected {expect}, got {actual}"
        )


def test_diatomic_curve_metrics(
    pred_ref_diatomic_curves: tuple[DiatomicCurves, DiatomicCurves],
) -> None:
    """Test full metrics calculation pipeline."""
    ref_curves, pred_curves = pred_ref_diatomic_curves
    ref_dists, pred_dists = ref_curves.distances, pred_curves.distances

    assert ref_dists == pytest.approx(pred_dists)

    # Test with default parameters (no force curves)
    metrics = diatomics.calc_diatomic_metrics(ref_curves, pred_curves)
    assert isinstance(metrics, dict)
    assert "H" in metrics
    metric_keys = [*metrics["H"]]
    assert set(metric_keys) >= {
        MbdKey.norm_auc,
        MbdKey.smoothness,
        MbdKey.tortuosity,
        MbdKey.energy_jump,
        MbdKey.energy_diff_flips,
        MbdKey.energy_grad_norm_max,
    }

    # Test with force curves
    metrics_with_forces = diatomics.calc_diatomic_metrics(ref_curves, pred_curves)
    assert isinstance(metrics_with_forces, dict)
    assert "H" in metrics_with_forces
    metric_keys_with_forces = [*metrics_with_forces["H"]]
    assert set(metric_keys_with_forces) >= {
        MbdKey.norm_auc,
        MbdKey.smoothness,
        MbdKey.tortuosity,
        MbdKey.conservation,
        MbdKey.energy_jump,
        MbdKey.energy_diff_flips,
        MbdKey.energy_grad_norm_max,
    }

    # Test with custom parameters
    custom_metrics: dict[str, dict[str, Any]] = {
        MbdKey.norm_auc: {"seps_range": (1.0, 4.0)},
    }
    custom_results = diatomics.calc_diatomic_metrics(
        ref_curves, pred_curves, metrics=custom_metrics
    )
    assert isinstance(custom_results, dict)
    assert "H" in custom_results
    for key in custom_metrics:
        assert key in custom_results["H"]
        assert custom_results["H"][key] != metrics["H"][key], f"{key=}"

    # Test with invalid metric name
    with pytest.raises(
        ValueError, match=re.escape("unknown_metrics={'invalid'}. Valid metrics=")
    ):
        diatomics.calc_diatomic_metrics(
            ref_curves, pred_curves, metrics={"invalid": {}}
        )


def test_write_metrics_to_yaml(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test writing diatomic metrics to YAML file."""
    # Create a temporary YAML file
    yaml_path = tmp_path / "model.yaml"
    yaml_path.write_text(text := "model: test_model")

    model = Model.mace_mp_0
    monkeypatch.setattr(Model, "yaml_path", yaml_path)

    # Test with empty metrics
    diatomics.write_metrics_to_yaml(model, {})
    assert yaml_path.read_text() == text

    # Test with valid metrics
    metrics: dict[str, dict[str, float]] = {
        "H": {
            MbdKey.smoothness: 1.0,
            MbdKey.tortuosity: 2.0,
            MbdKey.conservation: 3.0,
        },
        "He": {
            MbdKey.smoothness: 4.0,
            MbdKey.tortuosity: 5.0,
            MbdKey.conservation: 6.0,
        },
    }
    diatomics.write_metrics_to_yaml(model, metrics)

    # Check that metrics were written correctly
    yaml_content = yaml_path.read_text()
    assert "metrics:" in yaml_content
    assert "diatomics:" in yaml_content
    assert "smoothness: 2.5" in yaml_content  # mean of [1.0, 4.0]
    assert "tortuosity: 3.5" in yaml_content  # mean of [2.0, 5.0]
    assert "conservation: 4.5" in yaml_content  # mean of [3.0, 6.0]
