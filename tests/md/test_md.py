"""Tests for molecular dynamics rollout helpers and benchmark pipeline."""

import functools
import json
import os
import platform
import sys
from pathlib import Path

import ase.io
import h5py
import numpy as np
import pandas as pd
import pytest
from ase import Atoms
from ase.build import bulk
from ase.calculators.emt import EMT

from matbench_discovery import today
from matbench_discovery.md import (
    NpzValue,
    list_reference_systems,
    read_reference_trajectory,
    read_trajectory,
    run_md_benchmark,
    run_nvt_md,
    run_nvt_md_resumable,
    validate_pred_trajectory,
    write_reference_h5,
)
from matbench_discovery.metrics import md as md_metrics
from matbench_discovery.trajectory import Trajectory


def _make_reference_h5(tmp_path: Path, *, system: str = "bulkCu_300K_test") -> str:
    """Create a tiny single-system reference HDF5 with dt=0.25 fs and T=300 K."""
    atoms = bulk("Cu", cubic=True) * (2, 2, 2)
    reference = run_nvt_md(
        atoms, EMT(), temperature_kelvin=300, n_steps=80, record_interval=1
    )
    ref_file = str(tmp_path / "ref.h5")
    # dt = saved-frame cadence = time_step_fs * record_interval = 0.25 * 1
    write_reference_h5(
        ref_file, [(system, Trajectory.from_ase(reference), 0.25, 300.0)]
    )
    return ref_file


def _make_private_reference_h5(
    tmp_path: Path, *, system: str = "bulkCu_300K_test", atoms: Atoms | None = None
) -> str:
    """Create a tiny private reference HDF5 retaining energy and force labels."""
    atoms = bulk("Cu", cubic=True) * (2, 2, 2) if atoms is None else atoms
    reference = Trajectory.from_ase(
        run_nvt_md(atoms, EMT(), temperature_kelvin=300, n_steps=80, record_interval=1)
    )
    ref_file = str(tmp_path / "private-ref.h5")
    with h5py.File(ref_file, "w") as file:
        group = file.create_group(system)
        reference.write_to_h5_group(group)
        group.attrs["dt_fs"] = 0.25
        group.attrs["temperature_kelvin"] = 300.0
    return ref_file


def _run_tiny_benchmark(
    ref_file: str,
    out_dir: Path,
    *,
    private_ref_file: str | None = None,
    systems: list[str] | None = None,
    dry_run: bool = False,
) -> pd.DataFrame:
    """Run the tiny EMT benchmark configuration shared by pipeline tests."""
    return run_md_benchmark(
        calculator=EMT(),
        model_key="emt",
        out_dir=str(out_dir),
        ref_file=ref_file,
        private_ref_file=private_ref_file,
        systems=systems,
        n_steps=20,  # 20 // 1 + 1 = 21 frames (odd -> final skips checkpoint)
        time_step_fs=1,
        record_interval=1,
        dry_run=dry_run,
    )


def _fail_save_checkpoint_after(
    monkeypatch: pytest.MonkeyPatch, *, kills: int = 3
) -> dict[str, int]:
    """Simulate requeues after the configured number of checkpoint writes."""
    import matbench_discovery.md as md_mod

    original_savez = md_mod._atomic_savez  # noqa: SLF001
    save_count = {"count": 0}

    def killing_savez(path: str, **arrays: NpzValue) -> None:
        """Save arrays and interrupt the configured checkpoint writes."""
        original_savez(path, **arrays)
        if path.endswith(".ckpt.npz"):
            save_count["count"] += 1
            if save_count["count"] <= kills:
                raise RuntimeError("simulated timeout")

    monkeypatch.setattr(md_mod, "_atomic_savez", killing_savez)
    return save_count


def test_reference_h5_roundtrip(tmp_path: Path) -> None:
    """Reference HDF5 round-trips data and rejects missing or empty systems."""
    ref_file = str(tmp_path / "ref.h5")
    trajectory_a = Trajectory.from_ase(
        run_nvt_md(
            bulk("Cu", cubic=True),
            EMT(),
            temperature_kelvin=300,
            n_steps=6,
            record_interval=1,
        )
    )
    trajectory_b = Trajectory.from_ase(
        run_nvt_md(
            bulk("Al", cubic=True),
            EMT(),
            temperature_kelvin=500,
            n_steps=4,
            record_interval=1,
        )
    )
    write_reference_h5(
        ref_file,
        [
            ("sysA_300K", trajectory_a, 0.25, 300.0),
            ("sysB_500K", trajectory_b, 1.5, 500.0),
        ],
    )

    assert list_reference_systems(ref_file) == ["sysA_300K", "sysB_500K"]

    trajectory, dt_fs, temperature = read_reference_trajectory(ref_file, "sysB_500K")
    assert dt_fs == 1.5
    assert temperature == 500.0
    np.testing.assert_array_equal(trajectory.positions, trajectory_b.positions)
    assert trajectory.energy is None
    assert trajectory.forces is None
    np.testing.assert_array_equal(trajectory.stress, trajectory_b.stress)

    # max_frames only reads the leading frames (cheap head reads of huge references)
    capped, _dt_fs, _temperature = read_reference_trajectory(
        ref_file, "sysA_300K", max_frames=3
    )
    assert capped.n_frames == 3
    np.testing.assert_array_equal(capped.positions, trajectory_a.positions[:3])

    with pytest.raises(KeyError, match="not in reference"):
        read_reference_trajectory(ref_file, "does_not_exist")

    with pytest.raises(ValueError, match="no systems"):
        write_reference_h5(ref_file, [])


def test_run_nvt_md() -> None:
    """NVT rollout records evenly spaced frames with full results attached."""
    atoms = bulk("Cu", cubic=True) * (2, 2, 2)
    initial_positions = atoms.positions.copy()
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
    np.testing.assert_allclose(atoms.positions, initial_positions, rtol=0, atol=0)
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
    """Equal seeds reproduce trajectories while different seeds diverge."""
    atoms = bulk("Cu", cubic=True)

    frames_a = run_nvt_md(atoms, EMT(), temperature_kelvin=300, n_steps=10, seed=42)
    frames_b = run_nvt_md(atoms, EMT(), temperature_kelvin=300, n_steps=10, seed=42)
    frames_c = run_nvt_md(atoms, EMT(), temperature_kelvin=300, n_steps=10, seed=7)

    np.testing.assert_allclose(
        frames_a[-1].positions, frames_b[-1].positions, rtol=0, atol=0
    )
    assert np.abs(frames_a[-1].positions - frames_c[-1].positions).max() > 1e-6


def test_md_pipeline_end_to_end(tmp_path: Path) -> None:
    """Rollout, extxyz round-trip, and per-system evaluation compose."""
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
        n_rdf_bins=100,
    )

    # same model run at the same temperature: structure and dynamics should be
    # similar but not identical (different seeds), pressures close
    assert 0 < system_metrics["rdf_error"] < 50
    assert 0 < system_metrics["adf_error"] < 100
    assert 0 < system_metrics["vdos_error"] < 100
    assert 0 <= system_metrics["pressure_mae"] < 5
    assert 0 <= system_metrics["pressure_wasserstein"] < 5

    df_md = pd.DataFrame([system_metrics])
    model_metrics = md_metrics.calc_md_metrics(df_md)
    assert model_metrics["n_systems"] == 1
    assert "combined_score" not in model_metrics  # CMDS is site-computed, not stored


def test_validate_pred_trajectory(tmp_path: Path) -> None:
    """Complete rollouts validate while malformed variants return None."""
    atoms = bulk("Cu", cubic=True) * (2, 2, 2)
    frames = run_nvt_md(
        atoms, EMT(), temperature_kelvin=300, n_steps=40, record_interval=10
    )
    assert len(frames) == 5  # 40 // 10 + 1
    good = tmp_path / "pred.extxyz"
    ase.io.write(good, frames)

    def validate(
        path: Path, *, expected: int = 5, ref: Atoms = atoms
    ) -> list[Atoms] | None:
        """Validate a rollout against the expected frame and atom metadata."""
        return validate_pred_trajectory(
            str(path), expected_frames=expected, ref_atoms=ref, record_interval=10
        )

    valid = validate(good)
    assert valid is not None
    assert len(valid) == 5

    assert validate(good, expected=6) is None  # frame-count mismatch
    assert validate(tmp_path / "missing.extxyz") is None  # unreadable

    truncated = tmp_path / "trunc.extxyz"
    ase.io.write(truncated, frames[:3])
    assert validate(truncated) is None  # too few frames

    assert validate(good, ref=bulk("Al", cubic=True) * (2, 2, 2)) is None  # symbols

    for name, mutate in {
        "nan": lambda frames: frames[2].positions.__setitem__((0, 0), np.nan),
        "step": lambda frames: frames[2].info.__setitem__("md_step", 999),
    }.items():
        bad_frames = [frame.copy() for frame in frames]
        mutate(bad_frames)
        bad_path = tmp_path / f"{name}.extxyz"
        ase.io.write(bad_path, bad_frames)
        assert validate(bad_path) is None, f"{name} should fail validation"


def test_nvt_checkpoint_state_attrs() -> None:
    """ASE NVT must retain the private state attributes used by checkpoints."""
    from matbench_discovery.md import _init_nvt_dynamics, _nvt_state_attrs

    atoms = bulk("Cu", cubic=True) * (2, 2, 2)
    atoms.calc = EMT()
    dynamics = _init_nvt_dynamics(
        atoms, temperature_kelvin=300, time_step_fs=0.25, thermostat_time_scale_fs=25
    )
    state = _nvt_state_attrs(dynamics)  # raises RuntimeError if ASE renamed attrs
    assert len(state) == 4
    assert all(isinstance(array, np.ndarray) for array in state)


def _run_short_resumable(atoms: Atoms, out_path: str) -> list[Atoms]:
    """Run the short deterministic resumable rollout used by resume tests."""
    return run_nvt_md_resumable(
        atoms,
        EMT(),
        out_path=out_path,
        temperature_kelvin=300,
        n_steps=60,  # 60 // 10 + 1 = 7 frames (odd -> last frame skips checkpoint)
        record_interval=10,
        checkpoint_every_n_frames=2,
        progress_interval=0,
    )


def test_run_nvt_md_resume_reproduces(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Repeated resume reproduces an uninterrupted run frame-for-frame."""
    atoms = bulk("Cu", cubic=True) * (2, 2, 2)
    reference = _run_short_resumable(atoms, str(tmp_path / "ref.extxyz"))
    assert len(reference) == 7

    # simulate timeouts: raise after each of the first 3 checkpoint writes, then let
    # the run finish; the outer loop reruns (resuming) until it completes
    save_count = _fail_save_checkpoint_after(monkeypatch)

    out_path = str(tmp_path / "out.extxyz")
    resumed = None
    for _attempt in range(8):  # bounded resume loop
        try:
            resumed = _run_short_resumable(atoms, out_path)
            break
        except RuntimeError:
            continue
    assert resumed is not None, "resume loop did not complete"
    assert save_count["count"] >= 3, "test did not exercise interruptions"

    assert len(resumed) == len(reference)
    for frame, ref_frame in zip(resumed, reference, strict=True):
        assert frame.info["md_step"] == ref_frame.info["md_step"]
        np.testing.assert_allclose(
            frame.positions, ref_frame.positions, rtol=0, atol=1e-14
        )
        np.testing.assert_allclose(
            frame.get_velocities(), ref_frame.get_velocities(), rtol=0, atol=1e-14
        )
    # staged artifacts cleaned up on completion
    assert not os.path.isfile(f"{tmp_path}/out.part.extxyz")
    assert not os.path.isfile(f"{out_path}.ckpt.npz")


def test_run_nvt_md_resumable_ignores_corrupt_checkpoint(tmp_path: Path) -> None:
    """Corrupt checkpoints fail closed and trigger a clean rollout."""
    atoms = bulk("Cu", cubic=True) * (2, 2, 2)
    out_path = str(tmp_path / "out.extxyz")
    (tmp_path / "out.extxyz.ckpt.npz").write_bytes(b"not a real npz checkpoint")
    (tmp_path / "out.part.extxyz").write_text("garbage not extxyz\n")

    frames = run_nvt_md_resumable(
        atoms,
        EMT(),
        out_path=out_path,
        temperature_kelvin=300,
        n_steps=30,
        record_interval=10,
        progress_interval=0,
    )
    assert len(frames) == 4  # 30 // 10 + 1, recomputed cleanly
    assert os.path.isfile(out_path)


@pytest.mark.parametrize("dry_run", [True, False])
def test_run_md_benchmark(tmp_path: Path, *, dry_run: bool) -> None:
    """Dry and full runs return metrics while only full runs persist outputs."""
    ref_file = _make_reference_h5(tmp_path)
    out_dir = tmp_path / "out"

    def run(*, dry: bool = False) -> pd.DataFrame:
        """Run the benchmark in dry or persistent mode."""
        return _run_tiny_benchmark(ref_file, out_dir, dry_run=dry)

    df_md = run(dry=dry_run)
    assert list(df_md.index) == ["bulkCu_300K_test"]
    metric_columns = {
        "rdf_error",
        "adf_error",
        "vdos_error",
        "pressure_mae",
        "pressure_wasserstein",
        "pressure_error",
    }
    assert metric_columns <= set(df_md.columns)

    if dry_run:
        assert not out_dir.exists()  # dry run must not write outputs
        assert "run_time_sec" not in df_md  # in-memory rollout leaves no sidecar
        return

    assert (out_dir / f"{today}-emt-md-metrics.csv.gz").is_file()
    rollout = out_dir / "bulkCu_300K_test-nvt-emt.extxyz"
    assert rollout.is_file()
    # completed rollout leaves a run-info sidecar whose cost lands in the CSV row
    run_info_file = out_dir / "bulkCu_300K_test-nvt-emt-run-info.json"
    assert run_info_file.is_file()
    run_info = json.loads(run_info_file.read_text())
    # audit-only provenance stays in the sidecar (not copied into CSV rows)
    assert {"completed_at", "hostname", "versions"} <= set(run_info)
    assert run_info["versions"]["python"] == platform.python_version()
    assert df_md.loc["bulkCu_300K_test", "run_time_sec"] > 0
    assert isinstance(df_md.loc["bulkCu_300K_test", "hardware"], str)
    assert df_md.loc["bulkCu_300K_test", "n_atoms"] == 32  # 2x2x2 cubic Cu cell
    if sys.platform != "win32":  # Windows lacks the resource module
        assert df_md.loc["bulkCu_300K_test", "max_rss_gb"] > 0
    # second run reuses the existing rollout but still reports its recorded cost
    rollout_mtime = rollout.stat().st_mtime
    df_rerun = run()
    assert rollout.stat().st_mtime == rollout_mtime
    assert df_rerun.loc["bulkCu_300K_test", "run_time_sec"] == pytest.approx(
        df_md.loc["bulkCu_300K_test", "run_time_sec"]
    )


def test_run_md_benchmark_private_ref_adds_energy_force_rmse(tmp_path: Path) -> None:
    """Private labeled references enable energy and force diagnostics."""
    ref_file = _make_reference_h5(tmp_path)
    private_ref_file = _make_private_reference_h5(tmp_path)

    df_private_md = _run_tiny_benchmark(
        ref_file, tmp_path / "out", private_ref_file=private_ref_file, dry_run=True
    )

    assert {"energy_rmse", "force_rmse"} <= set(df_private_md.columns)
    assert np.isfinite(df_private_md.loc["bulkCu_300K_test", "energy_rmse"])
    assert np.isfinite(df_private_md.loc["bulkCu_300K_test", "force_rmse"])


@pytest.mark.parametrize(
    ("private_system", "private_atoms", "error_class", "match"),
    [
        # stale file: right group name but different (4 vs 32 atom) structure
        (
            "bulkCu_300K_test",
            bulk("Cu", cubic=True),
            ValueError,
            "different atoms than",
        ),
        # incomplete file: system missing entirely
        ("other_system", None, KeyError, "bulkCu_300K_test"),
    ],
)
def test_run_md_benchmark_rejects_bad_private_ref(
    tmp_path: Path,
    private_system: str,
    private_atoms: Atoms | None,
    error_class: type[Exception],
    match: str,
) -> None:
    """Stale or incomplete private references fail loudly."""
    ref_file = _make_reference_h5(tmp_path)
    private_ref_file = _make_private_reference_h5(
        tmp_path, system=private_system, atoms=private_atoms
    )

    with pytest.raises(error_class, match=match):
        _run_tiny_benchmark(
            ref_file,
            tmp_path / "out",
            private_ref_file=private_ref_file,
            dry_run=True,
        )


def test_run_md_benchmark_resume_matches_single_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Repeated benchmark resume matches one uninterrupted benchmark run."""
    import matbench_discovery.md as md_mod

    ref_file = _make_reference_h5(tmp_path)

    def run(out_dir: Path) -> pd.DataFrame:
        """Run the tiny benchmark with explicit args to keep ty types narrow."""
        return _run_tiny_benchmark(ref_file, out_dir)

    single = run(tmp_path / "single")

    # checkpoint every 2 frames so the short run actually checkpoints, and raise after
    # the first 3 checkpoint writes to mimic repeated job timeouts
    real_resumable = md_mod.run_nvt_md_resumable
    monkeypatch.setattr(
        md_mod,
        "run_nvt_md_resumable",
        functools.partial(real_resumable, checkpoint_every_n_frames=2),
    )
    save_count = _fail_save_checkpoint_after(monkeypatch)

    out_dir = tmp_path / "resumed"
    resumed = None
    for _attempt in range(8):
        try:
            resumed = run(out_dir)
            break
        except RuntimeError:
            continue
    assert resumed is not None, "resume loop did not complete"
    assert save_count["count"] >= 3, "test did not exercise interruptions"

    # wall time and memory high-water marks are process state, not trajectory
    # observables - deterministic reproduction doesn't apply to them
    numeric_columns = single.select_dtypes("number").columns.drop(
        ["run_time_sec", "max_rss_gb", "max_gpu_mem_gb"], errors="ignore"
    )
    for column in numeric_columns:
        np.testing.assert_allclose(
            resumed[column].to_numpy(),
            single[column].to_numpy(),
            rtol=0,
            atol=1e-9,
            equal_nan=True,
            err_msg=f"metric {column} differs after resume",
        )
    # completed rollout promoted, staged artifacts removed
    rollout = out_dir / "bulkCu_300K_test-nvt-emt.extxyz"
    assert rollout.is_file()
    assert not (out_dir / "bulkCu_300K_test-nvt-emt.part.extxyz").is_file()
    assert not (out_dir / "bulkCu_300K_test-nvt-emt.extxyz.ckpt.npz").is_file()
    assert len(read_trajectory(str(rollout))) == 21
    # rollout wall time accumulated across the interrupted sessions
    assert resumed.loc["bulkCu_300K_test", "run_time_sec"] > 0


@pytest.mark.parametrize("systems", [[], ["does_not_exist"]])
def test_run_md_benchmark_rejects_empty_or_unknown_systems(
    tmp_path: Path, systems: list[str]
) -> None:
    """Empty or unmatched system selections fail clearly."""
    ref_file = _make_reference_h5(tmp_path)

    with pytest.raises(ValueError, match="No systems found"):
        _run_tiny_benchmark(ref_file, tmp_path / "out", systems=systems)


@pytest.mark.parametrize("dt_factor", [5.0, 0.2])
def test_run_md_benchmark_rejects_mislabeled_reference_dt(
    tmp_path: Path, dt_factor: float
) -> None:
    """Bad reference dt_fs attributes fail before rollout."""
    ref_file = _make_reference_h5(tmp_path)
    with h5py.File(ref_file, "r+") as file:
        file["bulkCu_300K_test"].attrs["dt_fs"] *= dt_factor

    with pytest.raises(ValueError, match="mislabeled or stale"):
        _run_tiny_benchmark(ref_file, tmp_path / "out", dry_run=True)
