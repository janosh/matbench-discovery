"""Tests for molecular dynamics rollout helpers."""

from pathlib import Path

import ase.io
import numpy as np
import pandas as pd
import pytest
from ase.build import bulk
from ase.calculators.emt import EMT

from matbench_discovery.md import (
    extract_temperature,
    load_frame_intervals,
    read_trajectory,
    resolve_frame_interval,
    run_nvt_md,
)
from matbench_discovery.metrics import md as md_metrics


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("bulkAu_1500K_Kapil", 1500),
        ("MAPbBr3_300K", 300),
        ("sample_295.5K_xyz", 295.5),
        ("TiSe2_400K", 400),
        ("B4K2_600K", 600),  # K inside formula must not match
        ("no-temperature", None),
        ("400Kelvin", None),  # K followed by letters must not match
    ],
)
def test_extract_temperature(name: str, expected: float | None) -> None:
    """Temperatures should be parsed from system names with a <number>K token."""
    assert extract_temperature(name) == expected


def test_load_frame_intervals(tmp_path: Path) -> None:
    """Frame intervals should be keyed by '<system>_<temp>K'."""
    csv_path = tmp_path / "settings.csv"
    csv_path.write_text(
        "System,temperature,stride,dt,padding\n"
        "bulkAg,600,5,1,32\n"
        "anthracene,293,1,0.5,320\n"
        "bulkAg_600K_extra,600,1,9,32\n"
    )
    intervals = load_frame_intervals(str(csv_path))
    assert intervals == {
        "bulkAg_600K": 1.0,
        "anthracene_293K": 0.5,
        "bulkAg_600K_extra_600K": 9.0,
    }
    assert resolve_frame_interval("bulkAg_600K_Kapil", intervals) == 1.0
    assert resolve_frame_interval("bulkAg_600K_extra_600K_Kapil", intervals) == 9.0
    assert resolve_frame_interval("bulkAg_600K2_Kapil", intervals) is None


def test_run_nvt_md() -> None:
    """NVT rollout should record evenly spaced frames with full results attached."""
    atoms = bulk("Cu", cubic=True) * (2, 2, 2)
    init_positions = atoms.positions.copy()
    n_steps, record_interval = 50, 10

    frames = run_nvt_md(
        atoms,
        EMT(),
        temperature_kelvin=300,
        n_steps=n_steps,
        time_step_fs=1,
        record_interval=record_interval,
        seed=0,
    )

    assert len(frames) == n_steps // record_interval + 1
    assert [frame.info["md_step"] for frame in frames] == [0, 10, 20, 30, 40, 50]
    # input structure must not be modified
    np.testing.assert_allclose(atoms.positions, init_positions, rtol=0, atol=0)
    for frame in frames:
        assert np.isfinite(frame.get_potential_energy())
        assert frame.get_forces().shape == (len(atoms), 3)
        assert frame.get_stress(voigt=True).shape == (6,)
        # NVT: cell must stay fixed
        np.testing.assert_allclose(frame.cell.array, atoms.cell.array, rtol=0, atol=0)
    # velocities should be thermalized, i.e. non-zero
    assert np.abs(frames[0].get_velocities()).max() > 0
    # atoms should have moved between frames
    assert np.abs(frames[-1].positions - frames[0].positions).max() > 1e-3


def test_run_nvt_md_seed_reproducibility() -> None:
    """Same seed should give identical trajectories, different seeds different ones."""
    atoms = bulk("Cu", cubic=True)
    kwargs = dict(temperature_kelvin=300, n_steps=10, record_interval=5)

    frames_a = run_nvt_md(atoms, EMT(), seed=42, **kwargs)
    frames_b = run_nvt_md(atoms, EMT(), seed=42, **kwargs)
    frames_c = run_nvt_md(atoms, EMT(), seed=7, **kwargs)

    np.testing.assert_allclose(
        frames_a[-1].positions, frames_b[-1].positions, rtol=0, atol=0
    )
    assert np.abs(frames_a[-1].positions - frames_c[-1].positions).max() > 1e-6


def test_md_pipeline_end_to_end(tmp_path: Path) -> None:
    """Rollout -> extxyz round-trip -> full per-system evaluation should compose."""
    atoms = bulk("Cu", cubic=True) * (2, 2, 2)
    ref_traj = run_nvt_md(
        atoms, EMT(), temperature_kelvin=300, n_steps=100, record_interval=5, seed=0
    )
    pred_traj = run_nvt_md(
        atoms, EMT(), temperature_kelvin=300, n_steps=200, record_interval=10, seed=1
    )

    # round-trip through extxyz to verify energies/forces/stress survive file IO
    traj_path = tmp_path / "pred.extxyz"
    ase.io.write(traj_path, pred_traj)
    pred_traj_io = read_trajectory(str(traj_path))

    system_metrics = md_metrics.evaluate_md_system(
        ref_traj,
        pred_traj_io,
        ref_time_step_fs=5 * 0.25,
        pred_time_step_fs=10 * 0.25,
        calculator=EMT(),
        n_rdf_bins=100,
    )

    # same model run at the same temperature: structure and dynamics should be
    # similar but not identical (different seeds), pressures close
    assert 0 < system_metrics["rdf_error"] < 50
    assert 0 < system_metrics["vdos_error"] < 100
    assert 0 <= system_metrics["pressure_mae"] < 5
    assert 0 <= system_metrics["pressure_wasserstein"] < 5
    # reference labels come from the same EMT calculator -> zero RMSE
    assert system_metrics["energy_rmse"] == pytest.approx(0, abs=1e-12)
    assert system_metrics["force_rmse"] == pytest.approx(0, abs=1e-12)

    df_md = pd.DataFrame([system_metrics])
    model_metrics = md_metrics.calc_md_metrics(df_md)
    assert model_metrics["n_systems"] == 1
    assert model_metrics["combined_error"] >= 0
