"""Test diatomic curve metrics calculation functions."""

import json
import re
from collections.abc import Mapping
from pathlib import Path

import numpy as np
import pytest

from matbench_discovery.enums import MbdKey, Model
from matbench_discovery.metrics import diatomics
from matbench_discovery.metrics.diatomics import DiatomicCurve, DiatomicCurves

np_rng = np.random.default_rng(seed=0)


@pytest.fixture
def diatomics_model(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[Model, Path]:
    """Model.mace_mp_0 with a fresh temp YAML patched in as its yaml_path."""
    yaml_path = tmp_path / "model.yaml"
    yaml_path.write_text("model: test_model")
    monkeypatch.setattr(Model, "yaml_path", yaml_path)
    return Model.mace_mp_0, yaml_path


def test_diatomic_classes() -> None:
    """Test DiatomicCurve and DiatomicCurves initialization and validation."""
    distances = np.logspace(1, -1, 40).tolist()
    energies = np_rng.random(len(distances)).tolist()
    forces = np_rng.random((len(distances), 2, 3)).tolist()

    curve = DiatomicCurve(distances=distances, energies=energies, forces=forces)
    assert {*map(type, (curve.distances, curve.energies, curve.forces))} == {np.ndarray}
    for actual, expected in [
        (curve.distances, distances),
        (curve.energies, energies),
        (curve.forces, forces),
    ]:
        np.testing.assert_array_equal(actual, np.asarray(expected))

    data = {
        "distances": distances,
        "homo-nuclear": {"H": {"energies": energies, "forces": forces}},
        "hetero-nuclear": {"H-He": {"energies": energies, "forces": forces}},
    }
    curves = DiatomicCurves.from_dict(data)
    assert "homo-nuclear" in data
    assert isinstance(curves.homo_nuclear["H"], DiatomicCurve)
    h_he_curve = curves.hetero_nuclear.get("H-He")
    assert isinstance(h_he_curve, DiatomicCurve)
    np.testing.assert_array_equal(curves.distances, distances)
    np.testing.assert_array_equal(curves.homo_nuclear["H"].energies, energies)
    np.testing.assert_array_equal(h_he_curve.energies, energies)


def test_load_dft_reference_curves(tmp_path: Path) -> None:
    """DFT reference loader converts bundled JSON schema to DiatomicCurves."""
    forces = [
        [[0.1, 0, 0], [-0.1, 0, 0]],
        [[0.0, 0, 0], [0.0, 0, 0]],
    ]
    ref_path = tmp_path / "diatomics-dft.json"
    ref_path.write_text(
        json.dumps(
            {
                "PBE": {
                    "H-H": {
                        "distances": [0.7, 1.0],
                        "energies": [0.2, 0.0],
                        "forces": forces,
                    }
                }
            }
        ),
        encoding="utf-8",
    )

    curves = diatomics.load_dft_reference_curves(ref_path=f"{ref_path}")

    assert list(curves.homo_nuclear) == ["H-H"]
    curve = curves.homo_nuclear["H-H"]
    np.testing.assert_array_equal(curve.distances, np.array([0.7, 1.0]))
    np.testing.assert_array_equal(curve.energies, np.array([0.2, 0.0]))
    np.testing.assert_array_equal(curve.forces, np.asarray(forces))


def base_curve(
    distances: np.ndarray, *, equilibrium_distance: float = 1.5
) -> np.ndarray:
    """Simple Morse potential."""
    return 5 * (1 - np.exp(-2 * (distances - equilibrium_distance))) ** 2 - 5


def h_curves(
    distances: np.ndarray, energies: np.ndarray, forces: np.ndarray
) -> DiatomicCurves:
    """Wrap one H-H curve in DiatomicCurves."""
    curve = DiatomicCurve(distances=distances, energies=energies, forces=forces)
    return DiatomicCurves(distances=distances, homo_nuclear={"H": curve})


@pytest.mark.parametrize(
    (
        "equilibrium_distance",
        "energy_shift",
        "pbe_energy_mae",
        "energy_jump",
    ),
    [
        (1.5, 0, 0, 1.22),
        (1.5, 1, 0, 1.22),
        (2.0, 0, 868.04, 3.35),
    ],
    ids=["no_mod", "vertical_shift", "morse"],
)
def test_curve_shifts(
    equilibrium_distance: float,
    energy_shift: float,
    pbe_energy_mae: float,
    energy_jump: float,
) -> None:
    """Test metrics for various curve modifications."""
    distances = np.logspace(1, -1, 40)
    ref_energies = base_curve(distances)
    pred_energies = (
        base_curve(distances, equilibrium_distance=equilibrium_distance) + energy_shift
    )
    dummy_forces = np.zeros((len(distances), 2, 3))
    dummy_forces[:, 0, 0] = -np.gradient(ref_energies, distances)
    dummy_forces[:, 1, 0] = np.gradient(ref_energies, distances)

    ref_curves = h_curves(distances, ref_energies, dummy_forces)
    pred_curves = h_curves(distances, pred_energies, dummy_forces)
    metrics_out = diatomics.calc_diatomic_metrics(ref_curves, pred_curves)

    expected_metrics = {
        MbdKey.pbe_energy_mae: pbe_energy_mae,
        MbdKey.tortuosity: 1,
        MbdKey.energy_jump: energy_jump,
    }
    for metric_key, expect in expected_metrics.items():
        actual = metrics_out["H"][metric_key]
        assert actual == pytest.approx(expect, abs=0.5), f"{metric_key=}"


def test_diatomic_curve_metrics(
    pred_ref_diatomic_curves: tuple[DiatomicCurves, DiatomicCurves],
) -> None:
    """Test full metrics calculation pipeline."""
    ref_curves, pred_curves = pred_ref_diatomic_curves
    ref_dists, pred_dists = ref_curves.distances, pred_curves.distances

    assert ref_dists == pytest.approx(pred_dists)

    # all energy + force metrics are computed (DiatomicCurve always carries forces)
    metrics = diatomics.calc_diatomic_metrics(ref_curves, pred_curves)
    assert "H" in metrics
    assert set(metrics["H"]) >= {
        MbdKey.tortuosity,
        MbdKey.energy_jump,
        MbdKey.energy_diff_flips,
        MbdKey.pbe_energy_mae,
    }

    # Test with custom parameters
    custom_metrics: dict[str, dict[str, object]] = {
        MbdKey.pbe_wall_dist_mae: {"thresholds_ev": (1.0,)},
    }
    custom_results = diatomics.calc_diatomic_metrics(
        ref_curves, pred_curves, metrics=custom_metrics
    )
    assert set(custom_results["H"]) == set(custom_metrics)
    assert custom_results["H"][MbdKey.pbe_wall_dist_mae] >= 0

    # Test with interpolation parameter
    # Create a copy of ref_curves with slightly different distances
    modified_dists = ref_dists.copy() * 1.001  # 0.1% difference
    modified_ref_curves = DiatomicCurves(
        distances=modified_dists,
        homo_nuclear={
            elem: DiatomicCurve(
                distances=modified_dists, energies=curve.energies, forces=curve.forces
            )
            for elem, curve in ref_curves.homo_nuclear.items()
        },
        hetero_nuclear=ref_curves.hetero_nuclear,
    )

    # This should raise an error without interpolation
    with pytest.raises(
        ValueError,
        match="Reference and predicted distances must be same when interpolate=False",
    ):
        diatomics.calc_diatomic_metrics(
            modified_ref_curves, pred_curves, interpolate=False
        )

    for interpolate in (True, 50):
        interp_results = diatomics.calc_diatomic_metrics(
            modified_ref_curves, pred_curves, interpolate=interpolate
        )
        assert "H" in interp_results

    # Test with invalid metric name
    with pytest.raises(
        ValueError, match=re.escape("unknown_metrics={'invalid'}. Valid metrics=")
    ):
        diatomics.calc_diatomic_metrics(
            ref_curves, pred_curves, metrics={"invalid": {}}
        )


def test_write_metrics_to_yaml(
    diatomics_model: tuple[Model, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Test writing diatomic metrics to YAML file."""
    from matbench_discovery import ROOT

    model, yaml_path = diatomics_model
    pred_file = "models/mace/mace-mp-0/2025-02-13-test-diatomics.json.gz"
    pred_file_url = "https://figshare.com/files/fake-url-00000"

    def set_diatomics_metrics(metrics: Mapping[str, object]) -> None:
        """Patch Model.metrics with a diatomics block."""
        monkeypatch.setattr(Model, "metrics", {"diatomics": dict(metrics)})

    set_diatomics_metrics({"pred_file": pred_file, "pred_file_url": pred_file_url})

    result = diatomics.write_metrics_to_yaml(model, {})
    assert result == {"pred_file": pred_file, "pred_file_url": pred_file_url}
    yaml_content = yaml_path.read_text()
    assert "diatomics:" in yaml_content
    assert f"pred_file_url: {pred_file_url}" in yaml_content

    set_diatomics_metrics({"pred_file": pred_file})
    assert diatomics.write_metrics_to_yaml(model, {}) == {"pred_file": pred_file}

    new_pred_file = "models/mace/mace-mp-0/2026-06-28-diatomics.json.gz"
    set_diatomics_metrics({})
    assert diatomics.write_metrics_to_yaml(model, {}, pred_file_path=new_pred_file) == {
        "pred_file": new_pred_file
    }
    assert diatomics.write_metrics_to_yaml(
        model, {}, pred_file_path=f"{ROOT}/{new_pred_file}"
    ) == {"pred_file": new_pred_file}
    with pytest.raises(ValueError, match="must be inside repo root"):
        diatomics.write_metrics_to_yaml(
            model, {}, pred_file_path=f"{yaml_path.parent}/outside.json.gz"
        )

    set_diatomics_metrics({"pred_file": pred_file, "pred_file_url": pred_file_url})
    metrics_by_element: dict[str, dict[str, float]] = {
        "H": {
            MbdKey.energy_jump: 1.0,
            MbdKey.tortuosity: 2.0,
        },
        "He": {
            MbdKey.energy_jump: 4.0,
            MbdKey.tortuosity: 5.0,
        },
    }
    expected_metrics = {
        MbdKey.energy_jump: 2.5,
        MbdKey.tortuosity: 3.5,
    }
    result = diatomics.write_metrics_to_yaml(model, metrics_by_element)

    yaml_content = yaml_path.read_text()
    assert "metrics:" in yaml_content
    assert "diatomics:" in yaml_content
    for metric_key, metric_value in expected_metrics.items():
        assert f"{metric_key}: {metric_value}" in yaml_content

    assert result == {
        "pred_file": pred_file,
        "pred_file_url": pred_file_url,
        **expected_metrics,
    }

    # pred_file_path overrides the existing pred_file; run_metadata is recorded ahead of
    # the metric values
    set_diatomics_metrics({"pred_file_url": pred_file_url})
    result = diatomics.write_metrics_to_yaml(
        model,
        metrics_by_element,
        pred_file_path=new_pred_file,
        run_metadata={
            "hardware": "NVIDIA H100 80GB HBM3",
            "run_time_sec": 120.0,
            "excluded_formulas": ["He-He"],
            "invalid_key": "ignored",
        },
    )
    assert result == {
        "pred_file_url": pred_file_url,
        "pred_file": new_pred_file,
        "hardware": "NVIDIA H100 80GB HBM3",
        "run_time_sec": 120.0,
        "excluded_formulas": ["He-He"],
        **expected_metrics,
    }
    yaml_content = yaml_path.read_text()
    assert yaml_content.index("hardware:") < yaml_content.index("energy_jump:")

    # a recompute without run_metadata preserves the existing run metadata
    set_diatomics_metrics(result)
    recomputed = diatomics.write_metrics_to_yaml(model, metrics_by_element)
    assert recomputed["hardware"] == "NVIDIA H100 80GB HBM3"
    assert recomputed["run_time_sec"] == 120.0


def test_eval_window(monkeypatch: pytest.MonkeyPatch) -> None:
    """_eval_window returns element-specific bounds and handles missing radius data."""
    from ase.data import atomic_numbers, covalent_radii, vdw_alvarez

    from matbench_discovery.metrics.diatomics import _eval_window

    atomic_num_h = atomic_numbers["H"]
    r_min, r_max = _eval_window("H-H", 6.0)
    assert r_min == pytest.approx(0.9 * covalent_radii[atomic_num_h])
    assert r_max == pytest.approx(min(3.1 * vdw_alvarez.vdw_radii[atomic_num_h], 6.0))
    # seps_max caps r_max (3.1 * vdW(Cu) > 6)
    assert _eval_window("Cu-Cu", 6.0)[1] == pytest.approx(6.0)
    assert _eval_window("H-H", 2.5)[1] == pytest.approx(2.5)  # capped by seps_max

    atomic_num_og = atomic_numbers["Og"]
    r_min, r_max = _eval_window("Og-Og", 6.0)
    assert r_min == pytest.approx(0.9 * covalent_radii[atomic_num_og])
    assert r_max == pytest.approx(6.0)

    monkeypatch.setattr(diatomics, "covalent_radii", np.ones(10))
    monkeypatch.setattr(diatomics.vdw_alvarez, "vdw_radii", np.ones(10))
    assert _eval_window("Og-Og", 6.0) == (0.0, 6.0)


def test_window_excludes_deep_overlap() -> None:
    """Metrics ignore the unphysical deep-overlap region below 0.9*r_cov."""
    dists = np.linspace(0.1, 6.0, 60)  # spans below H's window (~0.28 Å) and above
    energies = np.exp(-dists)  # smooth, small gradient in-window
    energies[dists < 0.2] = 1e6  # huge spike in the excluded deep-overlap region
    forces = np.zeros((len(dists), 2, 3))
    pred = h_curves(dists, energies, forces)
    metrics = diatomics.calc_diatomic_metrics(ref_curves=None, pred_curves=pred)
    # the 1e6 spike is below H's r_min so the in-window smooth curve has no jump
    assert metrics["H"][MbdKey.energy_jump] == pytest.approx(0)


def test_write_metrics_drops_deprecated_and_handles_nan(
    diatomics_model: tuple[Model, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Recompute drops deprecated keys, skips NaN elements, unions per-element keys."""
    model, yaml_path = diatomics_model
    # existing block carries a deprecated metric (smoothness) + a url
    monkeypatch.setattr(
        Model,
        "metrics",
        {
            "diatomics": {
                "pred_file_url": "https://figshare.com/files/x",
                "smoothness": 9.9,
            }
        },
    )

    # H has a finite tortuosity + an extra (ref-only) metric; He's tortuosity is NaN
    metrics: dict[str, dict[str, float]] = {
        "H": {MbdKey.tortuosity: 2.0, MbdKey.pbe_energy_mae: 1.0},
        "He": {MbdKey.tortuosity: float("nan")},
    }
    result = diatomics.write_metrics_to_yaml(model, metrics)

    assert "smoothness" not in result  # deprecated key fully dropped
    assert "smoothness" not in yaml_path.read_text()
    assert result["pred_file_url"] == "https://figshare.com/files/x"  # url preserved
    assert result["tortuosity"] == 2.0  # mean over the one finite value
    assert result["pbe_energy_mae"] == 1.0  # unioned key present only on H

    # a metric that is NaN for every element is omitted, not written as .nan
    all_nan = diatomics.write_metrics_to_yaml(
        model, {"H": {MbdKey.tortuosity: float("nan")}}
    )
    assert "tortuosity" not in all_nan
    assert "tortuosity" not in yaml_path.read_text()
