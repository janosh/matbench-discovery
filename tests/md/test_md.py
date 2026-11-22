"""Tests for molecular dynamics rollout helpers and benchmark pipeline."""

import dataclasses
import functools
import itertools
import json
import os
import platform
import sys
from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np
import pandas as pd
import pytest
from kemio import EMT, Structure
from kemio.io import Frame, read_extxyz, write_extxyz

from matbench_discovery import today
from matbench_discovery.calculators import load_calculator_with_runtime
from matbench_discovery.md import (
    CHECKPOINT_SCHEMA,
    THERMOSTAT_TIME_SCALE_FS,
    NpzValue,
    NvtParams,
    list_reference_systems,
    read_reference_trajectory,
    run_info_path,
    run_md_benchmark,
    run_nvt_md,
    run_nvt_md_resumable,
    validate_pred_trajectory,
    write_reference_h5,
)
from matbench_discovery.metrics import md as md_metrics
from matbench_discovery.trajectory import Trajectory

# per-process run provenance of an EMT "model" that took 1.5 s to load
EMT_RUNTIME = load_calculator_with_runtime("emt")[1] | {"load_time_sec": 1.5}


def bulk(symbol: str, *, cubic: bool = True) -> Structure:
    """Conventional fcc cell (4 atoms) of Cu or Al, like ase.build.bulk(cubic=True)."""
    assert cubic
    a = {"Cu": 3.61, "Al": 4.05}[symbol]
    frac = [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]
    return Structure(a * np.eye(3), [symbol] * 4, frac, fractional=True)


def _run_cu_reference_rollout(atoms: Structure) -> list[Frame]:
    """80-step 300 K EMT rollout saving every 0.25 fs step (81 frames)."""
    return run_nvt_md(
        atoms, EMT(), temperature_kelvin=300, n_steps=80, record_interval=1, seed=0
    )


@pytest.fixture(scope="module")
def cu_emt_frames() -> list[Frame]:
    """Shared 32-atom Cu EMT rollout (81 frames, dt=0.25 fs); tests must not mutate."""
    return _run_cu_reference_rollout(bulk("Cu", cubic=True).make_supercell([2, 2, 2]))


def _make_reference_h5(
    tmp_path: Path,
    frames: list[Frame],
    *,
    system: str = "bulkCu_300K_test",
    labeled: bool = False,
) -> str:
    """Write frames as a 300 K, dt=0.25 fs reference HDF5, optionally labeled."""
    reference = Trajectory.from_frames(frames)
    ref_file = str(tmp_path / ("private-ref.h5" if labeled else "ref.h5"))
    # dt = saved-frame cadence = time_step_fs * record_interval = 0.25 * 1
    if labeled:  # write_reference_h5 strips labels, so write the group directly
        with h5py.File(ref_file, "w") as file:
            group = file.create_group(system)
            reference.write_to_h5_group(group)
            group.attrs.update(dt_fs=0.25, temperature_kelvin=300.0)
    else:
        write_reference_h5(ref_file, [(system, reference, 0.25, 300.0)])
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
        model=EMT(),
        model_key="emt",
        out_dir=str(out_dir),
        ref_file=ref_file,
        private_ref_file=private_ref_file,
        systems=systems,
        n_steps=20,  # 20 // 1 + 1 = 21 frames (odd -> final skips checkpoint)
        time_step_fs=1,
        record_interval=1,
        dry_run=dry_run,
        runtime=EMT_RUNTIME,
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


def _resume_until_complete[T](run: Callable[[], T]) -> T:
    """Rerun an interrupted (simulated-timeout) job until it completes."""
    for _attempt in range(8):  # bounded resume loop
        try:
            return run()
        except RuntimeError:
            continue
    raise AssertionError("resume loop did not complete")


def test_reference_h5_roundtrip(tmp_path: Path) -> None:
    """Reference HDF5 round-trips data and rejects missing or empty systems."""
    ref_file = str(tmp_path / "ref.h5")
    trajectory_a = Trajectory.from_frames(
        run_nvt_md(
            bulk("Cu", cubic=True),
            EMT(),
            temperature_kelvin=300,
            n_steps=6,
            record_interval=1,
        )
    )
    trajectory_b = Trajectory.from_frames(
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
    atoms = bulk("Cu", cubic=True).make_supercell([2, 2, 2])
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
    assert [frame.structure.info["md_step"] for frame in frames] == [
        0,
        10,
        20,
        30,
        40,
        50,
    ]
    # input structure must not be modified
    np.testing.assert_allclose(atoms.positions, initial_positions, rtol=0, atol=0)
    for frame in frames:
        assert frame.energy is not None
        assert np.isfinite(frame.energy)
        assert frame.forces is not None
        assert frame.forces.shape == (len(atoms), 3)
        assert frame.stress is not None
        assert frame.stress.shape == (3, 3)
        # NVT: cell must stay fixed
        np.testing.assert_allclose(
            frame.structure.lattice, atoms.lattice, rtol=0, atol=0
        )
    # momenta should be thermalized, i.e. non-zero
    assert np.abs(frames[0].structure.site_data["momenta"]).max() > 0
    # atoms should have moved between frames
    displacement = frames[-1].structure.positions - frames[0].structure.positions
    assert np.abs(displacement).max() > 1e-3


class StressFailingEMT:
    """EMT whose stress head raises ``error`` (forces/energy still work)."""

    def __init__(self, error: Exception) -> None:
        """Store the error raised whenever stress is requested."""
        self.error = error

    def compute(
        self, structures: list[Structure], properties: list[str]
    ) -> list[dict[str, object]]:
        """Delegate to EMT unless stress is requested."""
        if "stress" in properties:
            raise self.error
        return EMT().compute(structures, properties)


@pytest.mark.parametrize(
    "error",
    [NotImplementedError("no stress head"), RuntimeError("CUDA out of memory")],
)
def test_run_nvt_md_stress_probe(error: Exception) -> None:
    """Only NotImplementedError means "no stress"; real model errors propagate."""
    atoms = bulk("Cu", cubic=True)

    def run() -> list:
        """Two NVT steps with the stress-failing model."""
        return run_nvt_md(
            atoms,
            StressFailingEMT(error),
            temperature_kelvin=300,
            n_steps=2,
            record_interval=1,
            seed=0,
        )

    if isinstance(error, NotImplementedError):
        frames = run()
        assert len(frames) == 3
        assert all(frame.stress is None for frame in frames)
    else:
        with pytest.raises(RuntimeError, match="CUDA out of memory"):
            run()


def test_run_nvt_md_seed_reproducibility() -> None:
    """Equal seeds reproduce trajectories while different seeds diverge."""
    atoms = bulk("Cu", cubic=True)

    frames_a = run_nvt_md(atoms, EMT(), temperature_kelvin=300, n_steps=10, seed=42)
    frames_b = run_nvt_md(atoms, EMT(), temperature_kelvin=300, n_steps=10, seed=42)
    frames_c = run_nvt_md(atoms, EMT(), temperature_kelvin=300, n_steps=10, seed=7)

    last_a, last_b = frames_a[-1].structure, frames_b[-1].structure
    np.testing.assert_allclose(last_a.positions, last_b.positions, rtol=0, atol=0)
    assert np.abs(last_a.positions - frames_c[-1].structure.positions).max() > 1e-6


def test_evaluate_md_system(tmp_path: Path, cu_emt_frames: list[Frame]) -> None:
    """Self and 2x-strided (dt doubled) copies score ~0; stress-less pred gives NaN."""
    # pred = leading half of the rollout re-saved at 2x stride (41 -> 21 frames, 0.5 fs
    # apart), round-tripped through extxyz so energies/forces/stress survive file IO
    pred_path = tmp_path / "pred.extxyz"
    write_extxyz(str(pred_path), cu_emt_frames[:41:2])
    strided_pred = read_extxyz(str(pred_path))
    ref_dt = 0.25  # 1 saved frame per 0.25 fs MD step

    def evaluate(pred: list[Frame], pred_dt: float) -> dict[str, float]:
        """Score pred against the full 81-frame reference rollout."""
        return md_metrics.evaluate_md_system(
            cu_emt_frames,
            pred,
            ref_time_step_fs=ref_dt,
            pred_time_step_fs=pred_dt,
            n_rdf_bins=100,
        )

    self_metrics = evaluate(cu_emt_frames, ref_dt)
    # energy/force RMSEs need private labels, so they come from run_md_benchmark only
    assert set(self_metrics) == set(md_metrics.PER_SYSTEM_METRIC_COLS) - {
        "energy_rmse",
        "force_rmse",
    }
    # identical inputs go through identical arithmetic -> exactly zero
    assert self_metrics == dict.fromkeys(self_metrics, 0.0)

    # the strided copy samples a subset of the same time window, so errors are only
    # sampling noise: measured rdf 0.085%, adf 0.030%, vdos 1.6%, pressure_mae 7e-4
    # GPa. A wrong pred dt (rdf 10%) or untruncated ref (rdf 15%) far exceed the
    # bounds. pressure_error (histogram overlap of 21 vs 41 values) is too noisy here.
    strided_metrics = evaluate(strided_pred, 2 * ref_dt)
    for key, max_err in {
        "rdf_error": 0.5,
        "adf_error": 0.2,
        "vdos_error": 3,
        "pressure_mae": 0.01,
        "pressure_wasserstein": 0.01,
    }.items():
        assert 0 < strided_metrics[key] < max_err, f"{key}={strided_metrics[key]}"

    # pred frames without stress yield NaN pressure metrics, not errors
    no_stress_pred = [dataclasses.replace(frame, stress=None) for frame in strided_pred]
    no_stress = evaluate(no_stress_pred, 2 * ref_dt)
    assert no_stress["rdf_error"] == strided_metrics["rdf_error"]
    for key in ("pressure_mae", "pressure_wasserstein", "pressure_error"):
        assert np.isnan(no_stress[key]), key

    model_metrics = md_metrics.calc_md_metrics(pd.DataFrame([strided_metrics]))
    assert model_metrics["n_systems"] == 1
    assert "combined_score" not in model_metrics  # CMDS is site-computed, not stored


def test_validate_pred_trajectory(tmp_path: Path, cu_emt_frames: list[Frame]) -> None:
    """Complete rollouts validate while malformed variants return None."""
    atoms = cu_emt_frames[0].structure
    # first 40 MD steps recorded every 10 steps: md_step 0, 10, ..., 40
    frames = cu_emt_frames[:41:10]
    assert len(frames) == 5
    good = tmp_path / "pred.extxyz"
    write_extxyz(str(good), frames)
    params = NvtParams(
        temperature_kelvin=300,
        time_step_fs=0.25,
        record_interval=10,
        thermostat_time_scale_fs=25,
        seed=0,
    )

    def validate(
        path: Path, *, expected: int = 5, ref: Structure = atoms
    ) -> list[Frame] | None:
        """Validate a rollout against the expected frame and atom metadata."""
        return validate_pred_trajectory(
            str(path), expected_frames=expected, ref_atoms=ref, params=params
        )

    valid = validate(good)  # no run-info sidecar (legacy rollout) -> trusted
    assert valid is not None
    assert len(valid) == 5

    # a sidecar recording other rollout parameters forces a recompute, matching ones
    # reuse, and a corrupt sidecar fails loudly instead of skipping the check
    sidecar = Path(run_info_path(str(good)))
    for recorded, expected_valid in (
        (dataclasses.replace(params, time_step_fs=0.5), False),
        (dataclasses.replace(params, seed=1), False),
        (params, True),
    ):
        sidecar.write_text(json.dumps({"nvt_params": dataclasses.asdict(recorded)}))
        assert (validate(good) is not None) is expected_valid, f"{recorded=}"
    sidecar.write_text("{truncated")
    with pytest.raises(ValueError, match="Corrupt MD run-info sidecar"):
        validate(good)
    sidecar.unlink()

    assert validate(good, expected=6) is None  # frame-count mismatch
    assert validate(tmp_path / "missing.extxyz") is None  # unreadable

    truncated = tmp_path / "trunc.extxyz"
    write_extxyz(str(truncated), frames[:3])
    assert validate(truncated) is None  # too few frames

    assert (
        validate(good, ref=bulk("Al", cubic=True).make_supercell([2, 2, 2])) is None
    )  # symbols

    def with_nan(frame: Frame) -> Frame:
        structure = frame.structure
        momenta = np.array(structure.site_data["momenta"], dtype=float)
        momenta[0, 0] = np.nan
        site_data = {**structure.site_data, "momenta": momenta}
        return dataclasses.replace(
            frame, structure=structure.replace(site_data=site_data)
        )

    def with_bad_step(frame: Frame) -> Frame:
        structure = frame.structure.copy()
        structure.info["md_step"] = 999
        return dataclasses.replace(frame, structure=structure)

    for name, mutate in {"nan": with_nan, "step": with_bad_step}.items():
        bad_frames = [
            mutate(frame) if idx == 2 else frame for idx, frame in enumerate(frames)
        ]
        bad_path = tmp_path / f"{name}.extxyz"
        write_extxyz(str(bad_path), bad_frames)
        assert validate(bad_path) is None, f"{name} should fail validation"

    # non-finite positions (edited in the file: a Structure cannot hold them)
    lines = good.read_text().splitlines()
    first_atom = lines[2].split()
    lines[2] = " ".join([first_atom[0], "nan", *first_atom[2:]])
    nan_positions = tmp_path / "nan-positions.extxyz"
    nan_positions.write_text("\n".join(lines) + "\n")
    assert validate(nan_positions) is None


def _run_short_resumable(
    atoms: Structure, out_path: str, runtime: dict[str, object] = EMT_RUNTIME
) -> list[Frame]:
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
        runtime=runtime,
    )


def test_run_nvt_md_resume_reproduces(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Repeated resume reproduces an uninterrupted run and accumulates its cost."""
    import matbench_discovery.md as md_mod

    # every perf_counter call advances 1 s, so steady-state run time is deterministic
    clock = SimpleNamespace(perf_counter=itertools.count().__next__)
    monkeypatch.setattr(md_mod, "time", clock)
    atoms = bulk("Cu", cubic=True).make_supercell([2, 2, 2])
    ref_path = str(tmp_path / "ref.extxyz")
    reference = _run_short_resumable(atoms, ref_path)
    assert len(reference) == 7

    # simulate timeouts: raise after each of the first 3 checkpoint writes, then let
    # the run finish; the outer loop reruns (resuming) until it completes
    save_count = _fail_save_checkpoint_after(monkeypatch)

    out_path = str(tmp_path / "out.extxyz")
    resumed = _resume_until_complete(lambda: _run_short_resumable(atoms, out_path))
    assert save_count["count"] >= 3, "test did not exercise interruptions"

    assert len(resumed) == len(reference)
    # restored positions pass through fractional storage (~1 ulp); measured drift
    # after 30 steps is ~1e-14, so 1e-12 still flags any real state-restore error
    for frame, ref_frame in zip(resumed, reference, strict=True):
        assert frame.structure.info["md_step"] == ref_frame.structure.info["md_step"]
        np.testing.assert_allclose(
            frame.structure.positions, ref_frame.structure.positions, rtol=0, atol=1e-12
        )
        np.testing.assert_allclose(
            np.asarray(frame.structure.site_data["momenta"], dtype=float),
            np.asarray(ref_frame.structure.site_data["momenta"], dtype=float),
            rtol=0,
            atol=1e-12,
        )
    # staged artifacts cleaned up on completion
    assert not os.path.isfile(f"{tmp_path}/out.part.extxyz")
    assert not os.path.isfile(f"{out_path}.ckpt.npz")
    # every session's model load and process count accumulate across resumes
    run_info = json.loads(Path(run_info_path(out_path)).read_text())
    n_sessions = 4  # 3 simulated timeouts + the completing session
    assert run_info["n_processes"] == n_sessions
    assert run_info["load_time_sec"] == pytest.approx(1.5 * n_sessions)
    assert run_info["cost_source"] == "runner"
    assert run_info["hardware"] == EMT_RUNTIME["hardware"]
    # run time sums over all 4 sessions: equal to the uninterrupted run's, while
    # the last session alone took only 1 tick
    ref_info = json.loads(Path(run_info_path(ref_path)).read_text())
    assert run_info["run_time_sec"] == ref_info["run_time_sec"] == 4


def test_run_nvt_md_resumable_ignores_other_kemio_build(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture
) -> None:
    """A checkpoint written by a different kemio build is discarded, not resumed."""
    import matbench_discovery.md as md_mod

    atoms = bulk("Cu", cubic=True).make_supercell([2, 2, 2])
    out_path = str(tmp_path / "out.extxyz")
    _fail_save_checkpoint_after(monkeypatch, kills=1)
    with pytest.raises(RuntimeError, match="simulated timeout"):
        _run_short_resumable(atoms, out_path)
    assert os.path.isfile(f"{out_path}.ckpt.npz")

    monkeypatch.setattr(md_mod, "BUILD_ID", "0.1.0+rebuilt")
    frames = _run_short_resumable(atoms, out_path)
    assert "schema/kemio build mismatch" in capsys.readouterr().out
    assert len(frames) == 7  # recomputed from step 0


def test_run_nvt_md_resumable_skips_truncated_tail(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture
) -> None:
    """A tail cut off mid-frame after the checkpoint still resumes the prefix."""
    atoms = bulk("Cu", cubic=True).make_supercell([2, 2, 2])
    reference = _run_short_resumable(atoms, str(tmp_path / "ref.extxyz"))
    out_path = str(tmp_path / "out.extxyz")
    _fail_save_checkpoint_after(monkeypatch, kills=1)
    with pytest.raises(RuntimeError, match="simulated timeout"):
        _run_short_resumable(atoms, out_path)
    part = tmp_path / "out.part.extxyz"
    part.write_text(part.read_text() + "32\nLattice='7 0 0 0 7 0 0 0 7'\nCu 0.0 0")

    frames = _run_short_resumable(atoms, out_path)
    assert "resuming" in capsys.readouterr().out
    assert [frame.structure.info["md_step"] for frame in frames] == [
        frame.structure.info["md_step"] for frame in reference
    ]


def test_run_nvt_md_resumable_ignores_corrupt_checkpoint(tmp_path: Path) -> None:
    """Corrupt checkpoints fail closed and trigger a clean rollout."""
    atoms = bulk("Cu", cubic=True).make_supercell([2, 2, 2])
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
        runtime=EMT_RUNTIME,
    )
    assert len(frames) == 4  # 30 // 10 + 1, recomputed cleanly
    assert os.path.isfile(out_path)


@pytest.mark.parametrize(
    "changed",
    [
        None,
        {"seed": 1},
        {"temperature_kelvin": 301.0},
        {"time_step_fs": 0.5},
        {"schema": CHECKPOINT_SCHEMA + 1},
    ],
    ids=["unchanged", "seed", "temperature", "time_step", "schema"],
)
def test_load_checkpoint_rejects_changed_run_parameters(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, changed: dict[str, float] | None
) -> None:
    """Checkpoints resume only under the same seed, T, dt, schema and hardware."""
    import matbench_discovery.md as md_mod

    atoms = bulk("Cu", cubic=True).make_supercell([2, 2, 2])
    out_path = str(tmp_path / "out.extxyz")
    # interrupt right after the first checkpoint write to leave a resumable state
    _fail_save_checkpoint_after(monkeypatch, kills=1)
    with pytest.raises(RuntimeError, match="simulated timeout"):
        _run_short_resumable(atoms, out_path)
    ckpt_path = f"{out_path}.ckpt.npz"
    assert os.path.isfile(ckpt_path)

    changed = dict(changed or {})
    if "schema" in changed:  # a newer code version bumped the checkpoint layout
        new_schema = changed.pop("schema")
        monkeypatch.setattr(md_mod, "CHECKPOINT_SCHEMA", new_schema)
    # params of _run_short_resumable (run_nvt_md_resumable defaults otherwise)
    params = NvtParams(
        temperature_kelvin=300,
        time_step_fs=0.25,
        record_interval=10,
        thermostat_time_scale_fs=THERMOSTAT_TIME_SCALE_FS,
        seed=0,
    )
    params = dataclasses.replace(params, **changed)
    loaded = md_mod._load_checkpoint(  # noqa: SLF001
        ckpt_path, atoms=atoms, params=params, n_steps=60
    )
    if changed or md_mod.CHECKPOINT_SCHEMA != CHECKPOINT_SCHEMA:
        assert loaded is None, f"stale checkpoint resumed after {changed=}"
    else:  # control: the untouched checkpoint must resume, else the cases are vacuous
        assert loaded is not None
        # first checkpoint holds frames at steps 0 and 10 (checkpoint every 2 frames)
        assert int(loaded[0]["nsteps"]) == 10
        # ... but not on other hardware, leaving it for a matching resume
        other_gpu = EMT_RUNTIME | {"hardware": "NVIDIA H200"}
        with pytest.raises(ValueError, match="cannot resume under different run"):
            _run_short_resumable(atoms, out_path, runtime=other_gpu)
        assert os.path.isfile(ckpt_path)


def test_run_md_benchmark_dry_run_skips_unsupported_systems(
    tmp_path: Path, cu_emt_frames: list[Frame]
) -> None:
    """EMT dry runs roll out the first reference system with elements EMT supports."""
    indium_frames = [
        Frame(
            Structure(
                frame.structure.lattice,
                ["In"] * len(frame.structure),
                frame.structure.positions,
            )
        )
        for frame in cu_emt_frames
    ]
    ref_file = str(tmp_path / "ref.h5")
    write_reference_h5(
        ref_file,
        [
            ("InP_test", Trajectory.from_frames(indium_frames), 0.25, 300.0),
            ("bulkCu_300K_test", Trajectory.from_frames(cu_emt_frames), 0.25, 300.0),
        ],
    )
    df_md = _run_tiny_benchmark(ref_file, tmp_path / "out", dry_run=True)
    assert list(df_md.index) == ["bulkCu_300K_test"]
    with pytest.raises(ValueError, match=r"No system in \['InP_test'\] has elements"):
        _run_tiny_benchmark(
            ref_file, tmp_path / "out", systems=["InP_test"], dry_run=True
        )


@pytest.mark.parametrize("dry_run", [True, False])
def test_run_md_benchmark(
    tmp_path: Path, cu_emt_frames: list[Frame], *, dry_run: bool
) -> None:
    """Dry and full runs return metrics while only full runs persist outputs."""
    ref_file = _make_reference_h5(tmp_path, cu_emt_frames)
    out_dir = tmp_path / "out"

    def run(*, dry: bool = False, systems: list[str] | None = None) -> pd.DataFrame:
        """Run the benchmark in dry or persistent mode."""
        return _run_tiny_benchmark(ref_file, out_dir, systems=systems, dry_run=dry)

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
    assert {"nvt_params", "completed_at", "hostname", "versions"} <= set(run_info)
    assert run_info["versions"]["python"] == platform.python_version()
    # recorded, non-negative: cost rounds to 0.01 s, so a 20-step EMT rollout may be 0
    assert df_md.loc["bulkCu_300K_test", "run_time_sec"] >= 0
    assert isinstance(df_md.loc["bulkCu_300K_test", "hardware"], str)
    row = df_md.loc["bulkCu_300K_test"]
    # acceleration "none": CSV encoding of an empty flag list
    cost = (row["cost_source"], row["load_time_sec"], row["acceleration"])
    assert cost == ("runner", 1.5, "none")
    assert df_md.loc["bulkCu_300K_test", "n_atoms"] == 32  # 2x2x2 cubic Cu cell
    if sys.platform != "win32":  # Windows lacks the resource module
        assert df_md.loc["bulkCu_300K_test", "max_rss_gb"] > 0
    # second run reuses the existing rollout but still reports its recorded cost; an
    # explicit (duplicated) selection names the CSV after the matched systems
    rollout_mtime = rollout.stat().st_mtime
    df_rerun = run(systems=["bulkCu_300K_test"] * 2)
    assert rollout.stat().st_mtime == rollout_mtime
    assert (out_dir / f"{today}-emt-md-metrics-bulkCu_300K_test.csv.gz").is_file()
    assert df_rerun.loc["bulkCu_300K_test", "run_time_sec"] == pytest.approx(
        df_md.loc["bulkCu_300K_test", "run_time_sec"]
    )


@pytest.mark.parametrize(
    ("private_system", "private_atoms", "error_class", "match"),
    [
        ("bulkCu_300K_test", None, None, None),
        # stale file: right group name but different (4 vs 32 atom) structure
        ("bulkCu_300K_test", bulk("Cu", cubic=True), ValueError, "different atoms"),
        # incomplete file: system missing entirely
        ("other_system", None, KeyError, "bulkCu_300K_test"),
    ],
    ids=["matching", "stale", "incomplete"],
)
def test_run_md_benchmark_private_ref(
    tmp_path: Path,
    cu_emt_frames: list[Frame],
    private_system: str,
    private_atoms: Structure | None,
    error_class: type[Exception] | None,
    match: str | None,
) -> None:
    """Matching private labels add energy/force RMSEs; stale or incomplete fail."""
    ref_file = _make_reference_h5(tmp_path, cu_emt_frames)
    private_frames = (
        cu_emt_frames
        if private_atoms is None
        else _run_cu_reference_rollout(private_atoms)
    )
    private_ref_file = _make_reference_h5(
        tmp_path, private_frames, system=private_system, labeled=True
    )
    run = functools.partial(
        _run_tiny_benchmark,
        ref_file,
        tmp_path / "out",
        private_ref_file=private_ref_file,
        dry_run=True,
    )
    if error_class is not None:
        with pytest.raises(error_class, match=match):
            run()
        return

    df_private_md = run()
    # EMT scored against its own EMT labels: identical single points -> zero RMSE
    for col in ("energy_rmse", "force_rmse"):
        assert df_private_md.loc["bulkCu_300K_test", col] == pytest.approx(0, abs=1e-9)


def test_run_md_benchmark_resume_matches_single_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cu_emt_frames: list[Frame]
) -> None:
    """Repeated benchmark resume matches one uninterrupted benchmark run."""
    import matbench_discovery.md as md_mod

    ref_file = _make_reference_h5(tmp_path, cu_emt_frames)

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
    resumed = _resume_until_complete(lambda: run(out_dir))
    assert save_count["count"] >= 3, "test did not exercise interruptions"

    # wall time and memory high-water marks are process state, not trajectory
    # observables - deterministic reproduction doesn't apply to them
    process_cost_cols = ["run_time_sec", "load_time_sec", "n_processes"]
    numeric_columns = single.select_dtypes("number").columns.drop(
        [*process_cost_cols, "max_rss_gb", "max_gpu_mem_gb"], errors="ignore"
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
    assert len(read_extxyz(str(rollout))) == 21
    # rollout wall time accumulated across the interrupted sessions (rounded to 0.01 s)
    assert resumed.loc["bulkCu_300K_test", "run_time_sec"] >= 0


@pytest.mark.parametrize(
    ("systems", "match"),
    [
        ([], "No systems found"),
        (["does_not_exist"], r"Unknown systems \['does_not_exist'\]"),
        # a partly matching selection must not silently drop the unknown name
        (
            ["bulkCu_300K_test", "typo_system"],
            r"Unknown systems \['typo_system'\].*available: \['bulkCu_300K_test'\]",
        ),
    ],
)
def test_run_md_benchmark_rejects_empty_or_unknown_systems(
    tmp_path: Path, cu_emt_frames: list[Frame], systems: list[str], match: str
) -> None:
    """Empty selections and any unknown system names fail clearly."""
    ref_file = _make_reference_h5(tmp_path, cu_emt_frames)

    with pytest.raises(ValueError, match=match):
        _run_tiny_benchmark(ref_file, tmp_path / "out", systems=systems)


@pytest.mark.parametrize("dt_factor", [5.0, 0.2])
def test_run_md_benchmark_rejects_mislabeled_reference_dt(
    tmp_path: Path, cu_emt_frames: list[Frame], dt_factor: float
) -> None:
    """Bad reference dt_fs attributes fail before rollout."""
    ref_file = _make_reference_h5(tmp_path, cu_emt_frames)
    with h5py.File(ref_file, "r+") as file:
        file["bulkCu_300K_test"].attrs["dt_fs"] *= dt_factor

    with pytest.raises(ValueError, match="mislabeled or stale"):
        _run_tiny_benchmark(ref_file, tmp_path / "out", dry_run=True)
