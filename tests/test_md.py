"""Tests for molecular dynamics rollout helpers."""

import hashlib
import os
import shlex
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import Any

import ase.io
import h5py
import numpy as np
import pandas as pd
import pytest
from ase import Atoms
from ase.build import bulk
from ase.calculators.emt import EMT
from conftest import import_repo_script

from matbench_discovery import today
from matbench_discovery.calculators import CALCULATORS, load_calculator
from matbench_discovery.enums import Model
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


def test_reference_h5_roundtrip(tmp_path: Path) -> None:
    """write_reference_h5 then list_reference_systems / read_reference_trajectory
    round-trip preserves public frame data, dt_fs and temperature; missing systems
    raise; max_frames caps reads; an empty entries dict is rejected.
    """
    ref_file = str(tmp_path / "ref.h5")
    traj_a = Trajectory.from_ase(
        run_nvt_md(
            bulk("Cu", cubic=True),
            EMT(),
            temperature_kelvin=300,
            n_steps=6,
            record_interval=1,
        )
    )
    traj_b = Trajectory.from_ase(
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
        [("sysA_300K", traj_a, 0.25, 300.0), ("sysB_500K", traj_b, 1.5, 500.0)],
    )

    assert list_reference_systems(ref_file) == ["sysA_300K", "sysB_500K"]

    traj, dt_fs, temperature = read_reference_trajectory(ref_file, "sysB_500K")
    assert dt_fs == 1.5
    assert temperature == 500.0
    np.testing.assert_array_equal(traj.positions, traj_b.positions)
    assert traj.energy is None
    assert traj.forces is None
    np.testing.assert_array_equal(traj.stress, traj_b.stress)

    # max_frames only reads the leading frames (cheap head reads of huge references)
    capped, _dt, _temp = read_reference_trajectory(ref_file, "sysA_300K", max_frames=3)
    assert capped.n_frames == 3
    np.testing.assert_array_equal(capped.positions, traj_a.positions[:3])

    with pytest.raises(KeyError, match="not in reference"):
        read_reference_trajectory(ref_file, "does_not_exist")

    with pytest.raises(ValueError, match="no systems"):
        write_reference_h5(ref_file, [])


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

    frames_a = run_nvt_md(atoms, EMT(), temperature_kelvin=300, n_steps=10, seed=42)
    frames_b = run_nvt_md(atoms, EMT(), temperature_kelvin=300, n_steps=10, seed=42)
    frames_c = run_nvt_md(atoms, EMT(), temperature_kelvin=300, n_steps=10, seed=7)

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


def make_reference_h5(tmp_path: Path, *, system: str = "bulkCu_300K_test") -> str:
    """Tiny single-system reference HDF5 (dt=0.25 fs, T=300 K); returns the path."""
    atoms = bulk("Cu", cubic=True) * (2, 2, 2)
    ref = run_nvt_md(
        atoms, EMT(), temperature_kelvin=300, n_steps=80, record_interval=1
    )
    ref_file = str(tmp_path / "ref.h5")
    # dt = saved-frame cadence = time_step_fs * record_interval = 0.25 * 1
    write_reference_h5(ref_file, [(system, Trajectory.from_ase(ref), 0.25, 300.0)])
    return ref_file


def make_private_reference_h5(
    tmp_path: Path, *, system: str = "bulkCu_300K_test", atoms: Atoms | None = None
) -> str:
    """Tiny private reference HDF5 retaining energy/force labels for diagnostics."""
    if atoms is None:
        atoms = bulk("Cu", cubic=True) * (2, 2, 2)
    ref = Trajectory.from_ase(
        run_nvt_md(atoms, EMT(), temperature_kelvin=300, n_steps=80, record_interval=1)
    )
    ref_file = str(tmp_path / "private-ref.h5")
    with h5py.File(ref_file, "w") as file:
        group = file.create_group(system)
        ref.write_to_h5_group(group)
        group.attrs["dt_fs"] = 0.25
        group.attrs["temperature_kelvin"] = 300.0
    return ref_file


@pytest.mark.parametrize("dry_run", [True, False])
def test_run_md_benchmark(tmp_path: Path, *, dry_run: bool) -> None:
    """Both modes return one row per system with all metric columns. A dry run writes
    nothing; a full run writes a per-system CSV + a rollout it reuses on re-run.
    """
    ref_file = make_reference_h5(tmp_path)
    out_dir = tmp_path / "out"

    def run(*, dry: bool = False) -> pd.DataFrame:
        return run_md_benchmark(
            calculator=EMT(),
            model_key="emt",
            out_dir=str(out_dir),
            ref_file=ref_file,
            n_steps=20,
            time_step_fs=1,
            record_interval=1,
            dry_run=dry,
        )

    df_md = run(dry=dry_run)
    assert list(df_md.index) == ["bulkCu_300K_test"]
    metric_cols = {
        "rdf_error",
        "adf_error",
        "vdos_error",
        "pressure_mae",
        "pressure_wasserstein",
        "pressure_error",
    }
    assert metric_cols <= set(df_md.columns)

    if dry_run:
        assert not out_dir.exists()  # dry run must not write outputs
        return

    assert (out_dir / f"{today}-emt-md-metrics.csv.gz").is_file()
    rollout = out_dir / "bulkCu_300K_test-nvt-emt.extxyz"
    assert rollout.is_file()
    # second run reuses the existing rollout
    rollout_mtime = rollout.stat().st_mtime
    run()
    assert rollout.stat().st_mtime == rollout_mtime


def test_run_md_benchmark_private_ref_adds_energy_force_rmse(tmp_path: Path) -> None:
    """Private labeled reference file enables E/F diagnostics in benchmark rows."""
    ref_file = make_reference_h5(tmp_path)
    private_ref_file = make_private_reference_h5(tmp_path)

    df_private_md = run_md_benchmark(
        calculator=EMT(),
        model_key="emt",
        out_dir=str(tmp_path / "out"),
        ref_file=ref_file,
        private_ref_file=private_ref_file,
        n_steps=20,
        time_step_fs=1,
        record_interval=1,
        dry_run=True,
    )

    assert {"energy_rmse", "force_rmse"} <= set(df_private_md.columns)
    assert np.isfinite(df_private_md.loc["bulkCu_300K_test", "energy_rmse"])
    assert np.isfinite(df_private_md.loc["bulkCu_300K_test", "force_rmse"])


@pytest.mark.parametrize(
    ("private_kwargs", "err_cls", "match"),
    [
        # stale file: right group name but different (4 vs 32 atom) structure
        ({"atoms": bulk("Cu", cubic=True)}, ValueError, "different atoms than the"),
        # incomplete file: system missing entirely
        ({"system": "other_system"}, KeyError, "bulkCu_300K_test"),
    ],
)
def test_run_md_benchmark_rejects_bad_private_ref(
    tmp_path: Path,
    private_kwargs: dict[str, Any],
    err_cls: type[Exception],
    match: str,
) -> None:
    """A stale or incomplete private reference must fail loud, not yield
    plausible-but-wrong RMSEs or partial E/F coverage.
    """
    ref_file = make_reference_h5(tmp_path)
    private_ref_file = make_private_reference_h5(tmp_path, **private_kwargs)

    with pytest.raises(err_cls, match=match):
        run_md_benchmark(
            calculator=EMT(),
            model_key="emt",
            out_dir=str(tmp_path / "out"),
            ref_file=ref_file,
            private_ref_file=private_ref_file,
            n_steps=20,
            time_step_fs=1,
            record_interval=1,
            dry_run=True,
        )


def test_run_md_benchmark_resume_matches_single_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """End-to-end resume (the cluster timeout/requeue cycle): a rollout interrupted at
    several checkpoints and resumed across reruns of run_md_benchmark yields the same
    rollout and the same per-system metrics as a single uninterrupted run.
    """
    import functools

    import matbench_discovery.md as md_mod

    ref_file = make_reference_h5(tmp_path)

    def run(out_dir: Path) -> pd.DataFrame:
        """Run the tiny benchmark with explicit args (keeps ty from widening types)."""
        return run_md_benchmark(
            calculator=EMT(),
            model_key="emt",
            out_dir=str(out_dir),
            ref_file=ref_file,
            n_steps=20,  # 20 // 1 + 1 = 21 frames (odd -> final skips checkpoint)
            time_step_fs=1,
            record_interval=1,
        )

    single = run(tmp_path / "single")

    # checkpoint every 2 frames so the short run actually checkpoints, and raise after
    # the first 3 checkpoint writes to mimic repeated job timeouts
    real_resumable = md_mod.run_nvt_md_resumable
    monkeypatch.setattr(
        md_mod,
        "run_nvt_md_resumable",
        functools.partial(real_resumable, checkpoint_every_n_frames=2),
    )
    saves = _fail_save_checkpoint_after(monkeypatch)

    out_dir = tmp_path / "resumed"
    resumed = None
    for _attempt in range(8):
        try:
            resumed = run(out_dir)
            break
        except RuntimeError:
            continue
    assert resumed is not None, "resume loop did not complete"
    assert saves["n"] >= 3, "test did not exercise interruptions"

    numeric_cols = single.select_dtypes("number").columns
    for col in numeric_cols:
        np.testing.assert_allclose(
            resumed[col].to_numpy(),
            single[col].to_numpy(),
            rtol=0,
            atol=1e-9,
            equal_nan=True,
            err_msg=f"metric {col} differs after resume",
        )
    # completed rollout promoted, staged artifacts removed
    rollout = out_dir / "bulkCu_300K_test-nvt-emt.extxyz"
    assert rollout.is_file()
    assert not (out_dir / "bulkCu_300K_test-nvt-emt.part.extxyz").is_file()
    assert not (out_dir / "bulkCu_300K_test-nvt-emt.extxyz.ckpt.npz").is_file()
    assert len(read_trajectory(str(rollout))) == 21


def test_md_convert_resolve_settings() -> None:
    """resolve_settings does suffix-aware longest-prefix dir-name to key matching."""
    packer = import_repo_script(
        "md_convert_references_to_hdf5", "scripts/md_convert_references_to_hdf5.py"
    )

    settings = {
        "bulkAg_600K": (1.0, 600.0),
        "bulkAg_600K_extra_600K": (9.0, 600.0),
        "bulkCuAu_500K": (2.0, 500.0),
    }
    assert packer.resolve_settings("bulkAg_600K_Kapil", settings) == (1.0, 600.0)
    assert packer.resolve_settings("bulkAg_600K_extra_600K_Kapil", settings) == (
        9.0,
        600.0,
    )
    assert packer.resolve_settings("bulkAg_600K2_Kapil", settings) is None
    # '-' delimiter after the key must still resolve (caught a prod bug)
    assert packer.resolve_settings("bulkCuAu_500K-Artrith_VASP", settings) == (
        2.0,
        500.0,
    )


def test_md_convert_load_settings(tmp_path: Path) -> None:
    """load_settings reads dt_fs = the saved-frame interval = the CSV's dt column.
    The reference extxyz stores every integration step, so the stride column is NOT
    folded in (verified by equipartition; dt * stride underestimates temperature by
    stride**2 and collapses the vDOS frequency axis).
    """
    packer = import_repo_script(
        "md_convert_references_to_hdf5", "scripts/md_convert_references_to_hdf5.py"
    )

    csv = tmp_path / "settings.csv"
    csv.write_text(
        "System,temperature,stride,dt\nTiSe2,400,5,1\nanthracene,293,1,0.5\n"
    )
    settings = packer.load_settings(str(csv))
    assert settings["TiSe2_400K"] == (1.0, 400.0)  # dt = 1 fs (stride 5 not folded in)
    assert settings["anthracene_293K"] == (0.5, 293.0)  # dt = 0.5 fs


def test_md_convert_packs_reference(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The packer turns raw per-system extxyz + settings CSV into reference HDF5."""
    packer = import_repo_script(
        "md_convert_references_to_hdf5", "scripts/md_convert_references_to_hdf5.py"
    )

    frames = run_nvt_md(
        bulk("Cu", cubic=True) * (2, 2, 2),
        EMT(),
        temperature_kelvin=300,
        n_steps=8,
        record_interval=1,
    )
    system = "bulkCu_300K-Artrith_VASP"
    sys_dir = tmp_path / "raw" / system
    sys_dir.mkdir(parents=True)
    ase.io.write(sys_dir / "traj.extxyz", frames)
    settings_csv = tmp_path / "settings.csv"
    settings_csv.write_text("System,temperature,stride,dt\nbulkCu,300,4,0.25\n")
    out_h5 = str(tmp_path / "ref.h5")

    argv = [
        "pack",
        *("--ref-dir", str(tmp_path / "raw")),
        *("--settings-csv", str(settings_csv)),
        *("--out", out_h5),
    ]
    monkeypatch.setattr(sys, "argv", argv)
    assert packer.main() == 0

    assert list_reference_systems(out_h5) == [system]
    traj, dt_fs, temperature = read_reference_trajectory(out_h5, system)
    assert dt_fs == 0.25  # dt = 0.25 fs (stride 4 not folded into saved-frame interval)
    assert temperature == 300.0
    assert traj.n_frames == len(frames)
    # extxyz round-trips positions to print precision
    np.testing.assert_allclose(
        traj.positions[0], frames[0].positions, rtol=0, atol=1e-6
    )


@pytest.mark.parametrize("systems", [[], ["does_not_exist"]])
def test_run_md_benchmark_rejects_empty_or_unknown_systems(
    tmp_path: Path, systems: list[str]
) -> None:
    """An empty or unmatched --systems selection fails clearly."""
    ref_file = make_reference_h5(tmp_path)

    with pytest.raises(ValueError, match="No systems found"):
        run_md_benchmark(
            calculator=EMT(),
            model_key="emt",
            out_dir=str(tmp_path / "out"),
            ref_file=ref_file,
            systems=systems,
        )


def test_validate_pred_trajectory(tmp_path: Path) -> None:
    """A complete consistent rollout validates; truncated/short/NaN/wrong-step/
    mismatched-atoms files return None so the caller recomputes instead of silently
    evaluating a corrupt or shorter trajectory.
    """
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
        "nan": lambda fs: fs[2].positions.__setitem__((0, 0), np.nan),
        "step": lambda fs: fs[2].info.__setitem__("md_step", 999),
    }.items():
        bad_frames = [frame.copy() for frame in frames]
        mutate(bad_frames)
        bad_path = tmp_path / f"{name}.extxyz"
        ase.io.write(bad_path, bad_frames)
        assert validate(bad_path) is None, f"{name} should fail validation"


def test_nvt_checkpoint_state_attrs() -> None:
    """Guard: ASE NoseHooverChainNVT still exposes the private state attributes the
    checkpoint relies on, so an ASE upgrade that renames them fails loudly here.
    """
    from matbench_discovery.md import _init_nvt_dynamics, _nvt_state_attrs

    atoms = bulk("Cu", cubic=True) * (2, 2, 2)
    atoms.calc = EMT()
    dynamics = _init_nvt_dynamics(
        atoms, temperature_kelvin=300, time_step_fs=0.25, thermostat_time_scale_fs=25
    )
    state = _nvt_state_attrs(dynamics)  # raises RuntimeError if ASE renamed attrs
    assert len(state) == 4
    assert all(isinstance(array, np.ndarray) for array in state)


def _fail_save_checkpoint_after(
    monkeypatch: pytest.MonkeyPatch, *, kills: int = 3
) -> dict[str, int]:
    """Raise after the first ``kills`` checkpoint writes to simulate requeues."""
    import matbench_discovery.md as md_mod

    orig_savez = md_mod._atomic_savez  # noqa: SLF001
    saves = {"n": 0}

    def killing_savez(path: str, **arrays: NpzValue) -> None:
        orig_savez(path, **arrays)
        if path.endswith(".ckpt.npz"):
            saves["n"] += 1
            if saves["n"] <= kills:
                raise RuntimeError("simulated timeout")

    monkeypatch.setattr(md_mod, "_atomic_savez", killing_savez)
    return saves


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
    """A rollout interrupted at several checkpoints and resumed reproduces an
    uninterrupted run frame-for-frame (atol 1e-14), with no duplicate or missing
    frames.
    """
    atoms = bulk("Cu", cubic=True) * (2, 2, 2)
    reference = _run_short_resumable(atoms, str(tmp_path / "ref.extxyz"))
    assert len(reference) == 7

    # simulate timeouts: raise after each of the first 3 checkpoint writes, then let
    # the run finish; the outer loop reruns (resuming) until it completes
    saves = _fail_save_checkpoint_after(monkeypatch)

    out_path = str(tmp_path / "out.extxyz")
    resumed = None
    for _attempt in range(8):  # bounded resume loop
        try:
            resumed = _run_short_resumable(atoms, out_path)
            break
        except RuntimeError:
            continue
    assert resumed is not None, "resume loop did not complete"
    assert saves["n"] >= 3, "test did not exercise interruptions"

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
    """A corrupt/foreign checkpoint must be ignored (fail closed) and the rollout
    recomputed from scratch rather than restored into a bad state.
    """
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


def test_load_calculator(monkeypatch: pytest.MonkeyPatch) -> None:
    """Emt loads with no extra deps; dtype is ignored by models that don't declare it;
    unknown keys raise a helpful error.
    """
    from matbench_discovery import calculators

    assert isinstance(load_calculator("emt"), EMT)
    assert isinstance(load_calculator("emt", dtype="float32"), EMT)  # dtype ignored

    seen_dtype = ""

    def dtype_aware(device: str, dtype: str = "float64") -> EMT:  # noqa: ARG001
        nonlocal seen_dtype
        seen_dtype = dtype
        return EMT()

    monkeypatch.setitem(
        calculators.CALCULATORS, "mace_mp_0", calculators.CalcSpec(dtype_aware)
    )
    load_calculator("mace-mp-0", device="cpu", dtype="float32")
    assert seen_dtype == "float32"

    with pytest.raises(ValueError, match="Unknown model"):
        load_calculator("does-not-exist")


def test_download_checkpoint_replaces_zero_byte_cache(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """download_checkpoint redownloads zero-byte cached checkpoint files."""
    from matbench_discovery import calculators
    from matbench_discovery.enums import Model
    from matbench_discovery.remote import fetch

    url = "https://example.com/model.ckpt"
    monkeypatch.setattr(calculators, "CHECKPOINT_DIR", f"{tmp_path}")
    from_ref_mock = classmethod(
        lambda _cls, _model_key: SimpleNamespace(metadata={"checkpoint_url": url})
    )
    monkeypatch.setattr(Model, "from_ref", from_ref_mock)

    download_file_mock = lambda dest, _url, **_kwargs: Path(dest).write_bytes(  # noqa: E731
        b"checkpoint"
    )
    monkeypatch.setattr(fetch, "download_file", download_file_mock)
    url_hash = hashlib.sha256(url.encode()).hexdigest()[:12]
    dest = tmp_path / f"fake-{url_hash}.ckpt"
    dest.write_bytes(b"")

    actual_path = os.path.normpath(calculators.download_checkpoint("fake"))
    expected_path = os.path.normpath(f"{dest}")
    assert actual_path == expected_path
    assert dest.read_bytes() == b"checkpoint"


def test_run_md_cli_write_yaml_skips_non_submission_model(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`--write-yaml` for debug model emt must skip YAML write gracefully."""
    run_md = import_repo_script("run_md", "models/run_md.py")

    # avoid the real calculator + rollout pipeline; return metrics calc_md_metrics reads
    captured_kwargs: dict[str, object] = {}

    def fake_run_md_benchmark(**kwargs: object) -> pd.DataFrame:
        """Capture runner kwargs while returning minimal aggregateable metrics."""
        captured_kwargs.update(kwargs)
        return pd.DataFrame({"rdf_error": [1.0], "vdos_error": [2.0]})

    monkeypatch.setattr(run_md, "load_calculator", lambda *_a, **_k: EMT())
    monkeypatch.setattr(run_md, "run_md_benchmark", fake_run_md_benchmark)
    monkeypatch.setenv("MBD_MD_PRIVATE_REF_FILE", "/private/ref.h5")
    monkeypatch.setattr(sys, "argv", ["run_md", "--model", "emt", "--write-yaml"])

    assert run_md.main() == 0  # would raise ValueError without the non-enum guard
    assert captured_kwargs["private_ref_file"] == "/private/ref.h5"


@pytest.mark.parametrize(
    ("argv", "expected_cmd"),
    [
        (
            ["run_md", "--model", "mace_mp_0", "--write-yaml", "--systems", "bulkAu"],
            None,
        ),
        (
            ["run_md", "--print-cmd", "--model", "mace-mp-0"],
            CALCULATORS["mace_mp_0"].uv_run_cmd(
                "models/run_md.py", "--model", "mace_mp_0"
            ),
        ),
    ],
)
def test_run_md_cli_validation(
    argv: list[str],
    expected_cmd: list[str] | None,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """CLI rejects partial YAML writes and prints full canonical uv commands."""
    run_md = import_repo_script("run_md", "models/run_md.py")
    monkeypatch.setattr(sys, "argv", argv)

    if expected_cmd is None:
        with pytest.raises(SystemExit):  # parser.error exits before rollout pipeline
            run_md.main()
    else:
        assert run_md.main() == 0
        assert shlex.split(capsys.readouterr().out) == expected_cmd


def patch_eval_md(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    systems: list[str],
    model: Model = Model.mace_mp_0,
) -> ModuleType:
    """Import scripts/evals/md.py and patch it to use a tiny test reference set."""
    eval_md = import_repo_script("eval_md", "scripts/evals/md.py")
    monkeypatch.setattr(eval_md, "ROOT", str(tmp_path))
    monkeypatch.setattr(eval_md.cli_args, "models", [model])
    monkeypatch.setattr(eval_md, "default_md_reference_path", lambda: "ref.h5")
    monkeypatch.setattr(eval_md, "list_reference_systems", lambda _path: systems)
    return eval_md


def test_md_evals_skips_incomplete_coverage(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """scripts/evals/md.py must not write model YAML from incomplete coverage (only one
    of the three reference systems has a per-system CSV).
    """
    model = Model.mace_mp_0
    eval_md = patch_eval_md(tmp_path, monkeypatch, systems=["sysA", "sysB", "sysC"])

    def _fail_write(*_a: object, **_k: object) -> None:
        raise AssertionError("must not write YAML on incomplete coverage")

    monkeypatch.setattr(eval_md.md_metrics, "write_metrics_to_yaml", _fail_write)

    arch_dir = os.path.dirname(model.rel_path)
    md_dir = tmp_path / "models" / arch_dir / "2026-06-14-md-nvt"
    md_dir.mkdir(parents=True)
    pd.DataFrame({"system": ["sysA"], "rdf_error": [1.0]}).to_csv(
        md_dir / f"{model.name}-md-metrics-sysA.csv.gz", index=False
    )

    assert eval_md.main() == 1  # missing sysB/sysC -> skip, exit 1, no YAML written


def test_md_evals_skips_incomplete_fallback_csv(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The md_path fallback (a submitted combined CSV, no per-system files) must apply
    the same coverage guard so an incomplete submission can't corrupt the leaderboard.
    """
    model = Model.mace_mp_0
    eval_md = patch_eval_md(tmp_path, monkeypatch, systems=["sysA", "sysB", "sysC"])

    incomplete_csv = tmp_path / "combined.csv.gz"  # only 1 of 3 reference systems
    pd.DataFrame({"system": ["sysA"], "rdf_error": [1.0]}).to_csv(
        incomplete_csv, index=False
    )
    monkeypatch.setattr(
        type(model), "md_path", property(lambda _self: str(incomplete_csv))
    )

    def _fail_write(*_a: object, **_k: object) -> None:
        raise AssertionError("must not write YAML from an incomplete fallback CSV")

    monkeypatch.setattr(eval_md.md_metrics, "write_metrics_to_yaml", _fail_write)

    assert eval_md.main() == 1  # missing sysB/sysC -> skip, exit 1, no YAML written


def test_md_evals_handles_placeholder_md_metrics(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A model with complete per-system coverage but a 'not available' placeholder for
    metrics.md (a string, not a dict) must aggregate without crashing on md_yaml.get().
    """
    model = Model.mace_mp_0
    monkeypatch.setattr(
        type(model), "metrics", property(lambda _self: {"md": "not available"})
    )
    eval_md = patch_eval_md(tmp_path, monkeypatch, systems=["sysA"])

    writes: list[dict[str, object]] = []
    monkeypatch.setattr(
        eval_md.md_metrics,
        "write_metrics_to_yaml",
        lambda *_a, **kwargs: writes.append(kwargs),
    )

    arch_dir = os.path.dirname(model.rel_path)
    md_dir = tmp_path / "models" / arch_dir / "2026-06-14-md-nvt"
    md_dir.mkdir(parents=True)
    pd.DataFrame({"system": ["sysA"], "rdf_error": [1.0], "vdos_error": [2.0]}).to_csv(
        md_dir / f"{model.name}-md-metrics-sysA.csv.gz", index=False
    )

    assert eval_md.main() == 0  # placeholder md must not raise AttributeError
    assert len(writes) == 1
    assert writes[0]["pred_file_url"] is None  # no url from a placeholder md section


@pytest.mark.parametrize("model_key", [key for key in CALCULATORS if key != "emt"])
def test_md_model_keys_are_valid_enum_members(model_key: str) -> None:
    """Registered non-debug models must map to a Model so metrics write to YAML."""
    assert Model.from_ref(model_key).name == model_key
