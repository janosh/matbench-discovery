"""Tests for molecular dynamics metric helpers."""

from pathlib import Path
from types import SimpleNamespace

import ase.io
import numpy as np
import pandas as pd
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.singlepoint import SinglePointCalculator

from matbench_discovery.metrics import md


class InfoCalculator(Calculator):
    """ASE calculator returning predictions stored on the Atoms object."""

    # ASE exposes supported calculator results as a mutable class-level list.
    implemented_properties: list[str] = ["energy", "forces"]  # noqa: RUF012

    def calculate(
        self,
        atoms: Atoms | None = None,
        properties: list[str] | None = None,
        system_changes: list[str] = all_changes,
    ) -> None:
        """Populate energy and force predictions from frame metadata."""
        super().calculate(atoms, properties, system_changes)
        if atoms is None:
            raise ValueError("atoms must be provided")
        self.results["energy"] = float(atoms.info["pred_energy"])
        self.results["forces"] = np.asarray(atoms.arrays["pred_forces"], dtype=float)


class StressCalculator(Calculator):
    """ASE calculator returning a stress tensor stored on the Atoms object."""

    implemented_properties: list[str] = ["stress"]  # noqa: RUF012

    def calculate(
        self,
        atoms: Atoms | None = None,
        properties: list[str] | None = None,
        system_changes: list[str] = all_changes,
    ) -> None:
        """Populate stress predictions from frame metadata."""
        super().calculate(atoms, properties, system_changes)
        if atoms is None:
            raise ValueError("atoms must be provided")
        self.results["stress"] = np.asarray(atoms.info["pred_stress"], dtype=float)


def make_frame(
    formula: str,
    ref_energy: float,
    pred_energy: float,
    ref_forces: np.ndarray,
    pred_forces: np.ndarray,
) -> Atoms:
    """Create an ASE Atoms frame with reference and predicted labels."""
    atoms = Atoms(formula, positions=np.zeros((len(ref_forces), 3)))
    atoms.calc = SinglePointCalculator(atoms, energy=ref_energy, forces=ref_forces)
    atoms.info["pred_energy"] = pred_energy
    atoms.arrays["pred_forces"] = pred_forces
    return atoms


def make_rdf_frame(distance: float) -> Atoms:
    """Create a tiny periodic two-atom frame for RDF tests."""
    return Atoms(
        "Si2",
        positions=[(0, 0, 0), (distance, 0, 0)],
        cell=[8, 8, 8],
        pbc=True,
    )


def make_vdos_frame(displacement: float) -> Atoms:
    """Create a two-atom frame with one moving Cartesian coordinate."""
    return Atoms(
        "Si2",
        positions=[(0, 0, 0), (2 + displacement, 0, 0)],
        cell=[8, 8, 8],
        pbc=False,
    )


def make_pressure_frame(
    pressure_gpa: float,
    *,
    predicted_pressure_gpa: float | None = None,
    position_x: float = 0,
) -> Atoms:
    """Create a periodic frame with stored and optionally predicted pressure."""
    stress_value = -pressure_gpa / md.EV_PER_A3_TO_GPA
    stress = np.diag([stress_value, stress_value, stress_value])
    atoms = Atoms("Si", positions=[(position_x, 0, 0)], cell=[8, 8, 8], pbc=True)
    atoms.calc = SinglePointCalculator(atoms, stress=stress)
    if predicted_pressure_gpa is not None:
        predicted_stress = -predicted_pressure_gpa / md.EV_PER_A3_TO_GPA
        atoms.info["pred_stress"] = np.diag(
            [predicted_stress, predicted_stress, predicted_stress]
        )
    return atoms


def test_evaluate_reference_frames_and_metrics() -> None:
    """Energy and force RMSEs should be computed from reference MD frames."""
    frames = [
        make_frame(
            "H2",
            ref_energy=-2.0,
            pred_energy=-1.8,
            ref_forces=np.zeros((2, 3)),
            pred_forces=np.full((2, 3), 0.5),
        ),
        make_frame(
            "He",
            ref_energy=-1.0,
            pred_energy=-1.2,
            ref_forces=np.ones((1, 3)),
            pred_forces=np.full((1, 3), 1.5),
        ),
    ]

    df_eval = md.evaluate_reference_frames(
        frames,
        calculator=InfoCalculator(),
        model_key="test-model",
        input_file="sample_300K/traj.extxyz",
        system_name="sample_300K",
    )
    metrics = md.calc_md_reference_metrics(df_eval)

    assert df_eval["status"].tolist() == ["completed", "completed"]
    assert df_eval["temperature_k"].tolist() == [300, 300]
    assert metrics["energy_rmse"] == pytest.approx(np.sqrt((0.1**2 + 0.2**2) / 2))
    assert metrics["force_rmse"] == pytest.approx(0.5)
    assert metrics["n_frames"] == 2
    assert metrics["n_energy_frames"] == 2
    assert metrics["n_force_frames"] == 2


def test_calc_md_reference_metrics_handles_missing_values() -> None:
    """Aggregation should ignore missing labels and count failed rows."""
    df_eval = pd.DataFrame(
        {
            "frame_idx": [0, 1, 2],
            "status": ["completed", "failed", "missing_reference"],
            "energy_error_per_atom": [0.1, np.nan, np.nan],
            "force_sse": [3.0, np.nan, np.nan],
            "n_force_components": [12, 0, 0],
        }
    )

    metrics = md.calc_md_reference_metrics(df_eval)

    assert metrics["energy_rmse"] == pytest.approx(0.1)
    assert metrics["force_rmse"] == pytest.approx(0.5)
    assert metrics["n_frames"] == 3
    assert metrics["n_failed"] == 1


def test_calc_md_reference_metrics_empty_dataframe() -> None:
    """Empty prediction files should produce NaN errors and zero counts."""
    metrics = md.calc_md_reference_metrics(pd.DataFrame())

    assert np.isnan(metrics["energy_rmse"])
    assert np.isnan(metrics["force_rmse"])
    assert metrics["n_frames"] == 0
    assert metrics["n_energy_frames"] == 0
    assert metrics["n_force_frames"] == 0


def test_rdf_error_identical_trajectories() -> None:
    """Identical reference and MLIP trajectories should have zero RDF error."""
    frames = [make_rdf_frame(2.0), make_rdf_frame(2.1)]

    rdf = md.compute_rdf(frames, nbins=20, prefer_mdtraj=False)

    assert md.rdf_error(rdf, rdf) == pytest.approx(0)


def test_evaluate_rdf_dataset_requires_mlip_trajectories(tmp_path: Path) -> None:
    """RDF metrics should compare reference trajectories with MLIP rollouts."""
    ref_dir = tmp_path / "ref"
    ref_system_dir = ref_dir / "bulkSi_300K_Test"
    ref_system_dir.mkdir(parents=True)
    ref_frames = [make_rdf_frame(2.0), make_rdf_frame(2.1)]
    ase.io.write(ref_system_dir / "traj.extxyz", ref_frames)

    mlip_dir = tmp_path / "mlip"
    mlip_system_dir = mlip_dir / "bulkSi_300K_Test"
    mlip_system_dir.mkdir(parents=True)
    mlip_frames = [
        make_rdf_frame(2.0),
        make_rdf_frame(2.0),
        make_rdf_frame(2.1),
        make_rdf_frame(2.1),
    ]
    ase.io.write(mlip_system_dir / "nvt_test-model.extxyz", mlip_frames)

    df_eval, metrics = md.evaluate_rdf_dataset(
        model_key="test-model",
        ref_dir=ref_dir,
        mlip_dir=mlip_dir,
        output_dir=tmp_path / "out",
        config=md.MdRdfConfig(
            mode="same",
            nbins=20,
            ref_frame_dt_fs=2.0,
            mlip_frame_dt_fs=1.0,
            write_rdf_curves=False,
        ),
    )

    assert df_eval["status"].tolist() == ["completed"]
    assert df_eval["n_ref_frames"].tolist() == [2]
    assert df_eval["n_mlip_frames"].tolist() == [4]
    assert df_eval["matched_time_fs"].tolist() == [4]
    assert metrics["rdf_error"] == pytest.approx(0)
    assert metrics["n_rdf_systems"] == 1


def test_calc_md_rdf_metrics_requires_error_column() -> None:
    """RDF aggregation should not accept legacy RDF similarity columns."""
    with pytest.raises(ValueError, match="rdf_error or RDF_Error"):
        md.calc_md_rdf_metrics(pd.DataFrame({"rdf_similarity": [0.9]}))


def test_find_mlip_trajectory_accepts_ipi_layout(tmp_path: Path) -> None:
    """i-PI molecular-crystal runs should be discovered by system/model layout."""
    ipi_traj = (
        tmp_path / "anthracene_293K_Sharma_S" / "mace-mp-0" / "simulation.pos_0.extxyz"
    )
    ipi_traj.parent.mkdir(parents=True)
    ipi_traj.write_text("", encoding="utf-8")

    assert (
        md.find_mlip_trajectory(
            mlip_dir=tmp_path,
            system_name="anthracene_293K_Sharma_S",
            model_key="mace-mp-0",
        )
        == ipi_traj
    )


def test_vdos_error_identical_spectra() -> None:
    """Identical VDOS spectra should have zero source-normalized error."""
    frames = [
        make_vdos_frame(displacement) for displacement in (0.0, 0.2, 0.0, -0.2, 0.0)
    ]

    vdos = md.compute_vdos_hann(frames, frame_dt_fs=2.0, padding=2)

    assert md.vdos_error(vdos, vdos) == pytest.approx(0)


def test_evaluate_vdos_dataset_uses_source_settings_and_matching(
    tmp_path: Path,
) -> None:
    """VDOS settings should stride and match trajectories like the final scripts."""
    ref_dir = tmp_path / "ref"
    system_name = "bulkSi_300K_Test"
    ref_system_dir = ref_dir / system_name
    ref_system_dir.mkdir(parents=True)
    ref_frames = [
        make_vdos_frame(displacement)
        for displacement in (0.0, 0.1, 0.2, 0.1, 0.0, -0.1, -0.2, -0.1)
    ]
    ase.io.write(ref_system_dir / "traj.extxyz", ref_frames)

    mlip_dir = tmp_path / "mlip"
    mlip_system_dir = mlip_dir / system_name
    mlip_system_dir.mkdir(parents=True)
    mlip_frames = [
        *ref_frames[::2],
        make_vdos_frame(0.0),
        make_vdos_frame(0.2),
    ]
    ase.io.write(mlip_system_dir / "nvt_test-model.extxyz", mlip_frames)

    ref_settings = tmp_path / "vdos_settings_ref.csv"
    ref_settings.write_text(
        "System,temperature,stride,dt,padding\nbulkSi,300,2,1,2\n",
        encoding="utf-8",
    )
    mlip_settings = tmp_path / "vdos_settings_mlip.csv"
    mlip_settings.write_text(
        "System,temperature,stride,dt,padding\nbulkSi,300,1,2,2\n",
        encoding="utf-8",
    )

    output_dir = tmp_path / "out"
    df_eval, metrics = md.evaluate_vdos_dataset(
        model_key="test-model",
        ref_dir=ref_dir,
        mlip_dir=mlip_dir,
        output_dir=output_dir,
        config=md.MdVdosConfig(
            mode="same",
            ref_settings_path=ref_settings,
            mlip_settings_path=mlip_settings,
        ),
    )

    assert df_eval["status"].tolist() == ["completed"]
    assert df_eval["n_ref_frames"].tolist() == [4]
    assert df_eval["n_mlip_frames"].tolist() == [4]
    assert df_eval["ref_frame_dt_fs"].tolist() == [2]
    assert df_eval["mlip_frame_dt_fs"].tolist() == [2]
    assert metrics["vdos_error"] == pytest.approx(0)
    assert metrics["n_vdos_systems"] == 1
    curve_dir = output_dir / "vdos_results_hann_same_simulation_length"
    assert (
        curve_dir
        / "reference_matched"
        / system_name
        / "traj_stride2_match_nvt_test-model_vdos_hann.dat"
    ).is_file()
    assert (curve_dir / "mlip" / system_name / "nvt_test-model_vdos_hann.dat").is_file()


def test_calc_md_vdos_metrics_requires_error_column() -> None:
    """VDOS aggregation should not accept legacy VDOS similarity columns."""
    metrics = md.calc_md_vdos_metrics(
        pd.DataFrame({"final_mean_vdos_error": [0.25], "status": ["completed"]})
    )

    assert metrics["vdos_error"] == pytest.approx(0.25)
    with pytest.raises(ValueError, match="vdos_error"):
        md.calc_md_vdos_metrics(pd.DataFrame({"vdos_similarity": [0.9]}))


def test_pressure_conversion_and_identical_histogram_error() -> None:
    """Pressure conventions and identical distributions should match source behavior."""
    stress_value = -2.0 / md.EV_PER_A3_TO_GPA
    stress = np.diag([stress_value, stress_value, -10 / md.EV_PER_A3_TO_GPA])

    pressure_2d = md.pressure_from_stress_gpa(stress, pressure_mode="2d")
    scores = md.pressure_histogram_error(
        np.array([0.0, 1.0, 0.0, 1.0]),
        np.array([0.0, 1.0, 0.0, 1.0]),
        bins=4,
    )

    assert pressure_2d == pytest.approx(2)
    assert scores["pressure_error_percent"] == pytest.approx(0)


def test_evaluate_pressure_dataset_matches_time_and_uses_calculator(
    tmp_path: Path,
) -> None:
    """Pressure metrics should match times and evaluate model stress per frame."""
    ref_dir = tmp_path / "ref"
    system_name = "bulkSi_300K_Test"
    ref_system_dir = ref_dir / system_name
    ref_system_dir.mkdir(parents=True)
    ase.io.write(
        ref_system_dir / "traj.extxyz",
        [make_pressure_frame(value) for value in (0.0, 1.0, 0.0, 1.0)],
    )

    mlip_dir = tmp_path / "mlip"
    mlip_system_dir = mlip_dir / system_name
    mlip_system_dir.mkdir(parents=True)
    ase.io.write(
        mlip_system_dir / "nvt_test-model.extxyz",
        [
            make_pressure_frame(0, predicted_pressure_gpa=value, position_x=idx / 10)
            for idx, value in enumerate((0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0))
        ],
    )

    output_dir = tmp_path / "out"
    df_eval, metrics = md.evaluate_pressure_dataset(
        model_key="test-model",
        calculator=StressCalculator(),
        ref_dir=ref_dir,
        mlip_dir=mlip_dir,
        output_dir=output_dir,
        config=md.MdPressureConfig(
            bins=4,
            ref_frame_dt_fs=2.0,
            mlip_frame_dt_fs=1.0,
        ),
    )

    assert df_eval["status"].tolist() == ["completed"]
    assert df_eval["n_ref_frames"].tolist() == [4]
    assert df_eval["n_mlip_frames"].tolist() == [8]
    assert df_eval["matched_time_fs"].tolist() == [8]
    assert metrics["pressure_error_percent"] == pytest.approx(0)
    assert metrics["n_pressure_systems"] == 1
    values_dir = output_dir / "pressure_results_same_simulation_length"
    assert (
        values_dir / "reference" / f"{system_name}__test-model_pressure_per_frame.csv"
    ).is_file()
    assert (
        values_dir / "mlip" / f"{system_name}__test-model_pressure_per_frame.csv"
    ).is_file()


def test_calc_md_pressure_metrics_requires_error_column() -> None:
    """Pressure aggregation should not accept legacy pressure similarity columns."""
    with pytest.raises(ValueError, match="pressure_error_percent"):
        md.calc_md_pressure_metrics(pd.DataFrame({"pressure_similarity": [0.9]}))


def test_calc_md_combined_error_uses_fraction_scale() -> None:
    """Combined error should average all components on a 0-to-1 scale."""
    metrics = {
        "rdf_error": 10.0,
        "vdos_error": 0.4,
        "pressure_error_percent": 25.0,
    }

    assert md.calc_md_combined_error(metrics) == pytest.approx(0.25)
    assert md.calc_md_combined_error({"rdf_error": 10.0}) is None
    assert md.calc_md_combined_error(metrics | {"vdos_error": float("nan")}) is None
    assert (
        md.calc_md_combined_error(metrics | {"pressure_error_percent": float("nan")})
        is None
    )


def test_md_prediction_bundle_recomputes_combined_metrics(tmp_path: Path) -> None:
    """A bundled MD prediction artifact should reproduce all reported scores."""
    pred_path = md.write_md_prediction_bundle(
        model_key="test-model",
        output_dir=tmp_path,
        metric_frames={
            "reference": pd.DataFrame(
                {
                    "energy_error_per_atom": [0.2],
                    "force_sse": [3.0],
                    "n_force_components": [12],
                    "status": ["completed"],
                }
            ),
            "rdf": pd.DataFrame({"rdf_error": [10.0], "status": ["completed"]}),
            "vdos": pd.DataFrame({"vdos_error": [0.4], "status": ["completed"]}),
            "pressure": pd.DataFrame(
                {"pressure_error_percent": [25.0], "status": ["completed"]}
            ),
        },
    )

    bundled_rows = pd.read_csv(pred_path)
    metrics = md.calc_md_metrics(bundled_rows)

    assert pred_path.name == "md_predictions_bundle_test-model.csv"
    assert set(bundled_rows["metric_kind"]) == {"reference", "rdf", "vdos", "pressure"}
    assert metrics["energy_rmse"] == pytest.approx(0.2)
    assert metrics["force_rmse"] == pytest.approx(0.5)
    assert metrics["rdf_error"] == pytest.approx(10.0)
    assert metrics["vdos_error"] == pytest.approx(0.4)
    assert metrics["pressure_error_percent"] == pytest.approx(25)
    assert metrics["combined_error"] == pytest.approx(0.25)


def test_write_md_metrics_to_yaml(tmp_path: Path) -> None:
    """MD metrics should be written under metrics.md."""
    yaml_path = tmp_path / "model.yml"
    yaml_path.write_text("metrics:\n  md: not available\n", encoding="utf-8")
    model = SimpleNamespace(yaml_path=str(yaml_path))

    md.write_metrics_to_yaml(
        model,  # ty: ignore[invalid-argument-type]
        {
            "energy_rmse": 0.123456,
            "force_rmse": 0.654321,
            "rdf_error": 1.2345,
            "vdos_error": 0.234568,
            "vdos_similarity": 0.765432,
            "vdos_similarity_error": 0.234568,
            "pressure_similarity": 0.9234567,
            "pressure_error_fraction": 0.0765433,
            "pressure_similarity_percent": 92.34567,
            "pressure_error_percent": 7.65433,
            "combined_similarity": 0.8922,
            "n_frames": 10,
            "n_energy_frames": 10,
            "n_force_frames": 9,
            "n_rdf_systems": 3,
            "n_vdos_systems": 3,
            "n_pressure_systems": 3,
            "n_failed": 1,
            "n_skipped": 2,
        },
        pred_file_path="models/test/md_reference.csv",
    )

    text = yaml_path.read_text(encoding="utf-8")
    assert "md:" in text
    assert "energy_rmse: 0.1235" in text
    assert "force_rmse: 0.6543" in text
    assert "rdf_error: 1.2345" in text
    assert "vdos_error: 0.2346" in text
    assert "vdos_similarity:" not in text
    assert "vdos_similarity_error:" not in text
    assert "pressure_error_percent: 7.6543" in text
    assert "pressure_similarity:" not in text
    assert "pressure_similarity_percent:" not in text
    assert "pressure_error_fraction:" not in text
    assert "combined_error: 0.1078" in text
    assert "combined_similarity:" not in text
    assert "pred_file: models/test/md_reference.csv" in text
