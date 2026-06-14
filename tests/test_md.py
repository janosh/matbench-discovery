"""Tests for molecular dynamics rollout helpers."""

import os
from pathlib import Path

import ase.io
import numpy as np
import pandas as pd
import pytest
from ase import Atoms
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.md.nose_hoover_chain import NoseHooverChainNVT

from matbench_discovery.enums import Model
from matbench_discovery.md import (
    extract_temperature,
    load_frame_intervals,
    read_trajectory,
    resolve_frame_interval,
    run_md_benchmark,
    run_nvt_md,
    validate_pred_trajectory,
)
from matbench_discovery.md_models import MD_MODELS, load_calculator
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
        "sample,295.5,1,0.75,32\n"
        "bulkAg_600K_extra,600,1,9,32\n"
        "bulkCuAu,500,5,2,32\n"
    )
    intervals = load_frame_intervals(str(csv_path))
    assert intervals == {
        "bulkAg_600K": 1.0,
        "anthracene_293K": 0.5,
        "sample_295.5K": 0.75,
        "bulkAg_600K_extra_600K": 9.0,
        "bulkCuAu_500K": 2.0,
    }
    assert resolve_frame_interval("bulkAg_600K_Kapil", intervals) == 1.0
    assert resolve_frame_interval("sample_295.5K_Kapil", intervals) == 0.75
    assert resolve_frame_interval("bulkAg_600K_extra_600K_Kapil", intervals) == 9.0
    assert resolve_frame_interval("bulkAg_600K2_Kapil", intervals) is None
    # dir names may use '-' (not '_') as the delimiter after the '<system>_<temp>K'
    # key, e.g. bulkCuAu_500K-Artrith_VASP - must still resolve (caught a prod bug)
    assert resolve_frame_interval("bulkCuAu_500K-Artrith_VASP", intervals) == 2.0


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


def make_reference_dir(tmp_path: Path) -> tuple[str, str]:
    """Write a tiny CFPMD-style reference dir + settings CSV, return (ref_dir, csv)."""
    atoms = bulk("Cu", cubic=True) * (2, 2, 2)
    ref = run_nvt_md(
        atoms, EMT(), temperature_kelvin=300, n_steps=80, record_interval=1
    )
    sys_dir = tmp_path / "ref" / "bulkCu_300K_test"
    sys_dir.mkdir(parents=True)
    ase.io.write(sys_dir / "traj.extxyz", ref)
    csv = tmp_path / "settings.csv"
    # dt = saved-frame cadence = time_step_fs * record_interval = 0.25 * 1
    csv.write_text("System,temperature,stride,dt,padding\nbulkCu,300,1,0.25,32\n")
    return str(tmp_path / "ref"), str(csv)


@pytest.mark.parametrize("dry_run", [True, False])
def test_run_md_benchmark(tmp_path: Path, *, dry_run: bool) -> None:
    """Both modes return one row per system with all metric columns. A dry run writes
    nothing; a full run writes a per-system CSV + a rollout it reuses on re-run.
    """
    ref_dir, settings_csv = make_reference_dir(tmp_path)
    out_dir = tmp_path / "out"

    def run(*, dry: bool = False) -> pd.DataFrame:
        return run_md_benchmark(
            calculator=EMT(),
            model_key="emt",
            out_dir=str(out_dir),
            ref_dir=ref_dir,
            settings_csv=settings_csv,
            n_steps=20,
            time_step_fs=1,
            record_interval=1,
            dry_run=dry,
        )

    df_md = run(dry=dry_run)
    assert list(df_md.index) == ["bulkCu_300K_test"]
    metric_cols = {"rdf_error", "vdos_error", "energy_rmse", "force_rmse"}
    assert metric_cols <= set(df_md.columns)

    if dry_run:
        assert not out_dir.exists()  # dry run must not write outputs
        return

    assert (out_dir / "emt-md-metrics.csv.gz").is_file()
    rollout = out_dir / "bulkCu_300K_test-nvt-emt.extxyz"
    assert rollout.is_file()
    # second run reuses the existing rollout instead of recomputing
    mtime = rollout.stat().st_mtime
    run()
    assert rollout.stat().st_mtime == mtime


def test_run_md_benchmark_resume_matches_single_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """End-to-end resume (the cluster timeout/requeue cycle): a rollout interrupted at
    several checkpoints and resumed across reruns of run_md_benchmark yields the same
    rollout and the same per-system metrics as a single uninterrupted run.
    """
    import functools

    import matbench_discovery.md as md_mod

    ref_dir, settings_csv = make_reference_dir(tmp_path)

    def run(out_dir: Path) -> pd.DataFrame:
        """Run the tiny benchmark with explicit args (keeps ty from widening types)."""
        return run_md_benchmark(
            calculator=EMT(),
            model_key="emt",
            out_dir=str(out_dir),
            ref_dir=ref_dir,
            settings_csv=settings_csv,
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


def test_run_md_benchmark_rejects_invalid_inputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Invalid benchmark inputs should fail clearly before downstream errors.

    Covers an empty system selection (would otherwise become an empty-DataFrame
    KeyError) and an empty reference trajectory (would otherwise index frame 0).
    """
    import matbench_discovery.md as md_mod

    ref_dir, settings_csv = make_reference_dir(tmp_path)

    with pytest.raises(ValueError, match="No system directories found"):
        run_md_benchmark(
            calculator=EMT(),
            model_key="emt",
            out_dir=str(tmp_path / "out"),
            ref_dir=ref_dir,
            settings_csv=settings_csv,
            systems=["does_not_exist"],
        )

    monkeypatch.setattr(md_mod, "read_trajectory", lambda *_args, **_kwargs: [])
    with pytest.raises(ValueError, match="empty or corrupted"):
        run_md_benchmark(
            calculator=EMT(),
            model_key="emt",
            out_dir=str(tmp_path / "out"),
            ref_dir=ref_dir,
            settings_csv=settings_csv,
            n_steps=20,
            time_step_fs=1,
            record_interval=1,
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
    """Make _save_checkpoint raise after its first ``kills`` writes (simulating cluster
    timeouts that requeue the job), then succeed. Returns a dict tracking write count.
    """
    import matbench_discovery.md as md_mod

    orig_save = md_mod._save_checkpoint  # noqa: SLF001
    saves = {"n": 0}

    def killing_save(
        ckpt_path: str,
        dynamics: NoseHooverChainNVT,
        *,
        n_steps: int,
        n_frames: int,
        record_interval: int,
        time_step_fs: float,
        temperature_kelvin: float,
        thermostat_time_scale_fs: float,
        seed: int,
    ) -> None:
        orig_save(
            ckpt_path,
            dynamics,
            n_steps=n_steps,
            n_frames=n_frames,
            record_interval=record_interval,
            time_step_fs=time_step_fs,
            temperature_kelvin=temperature_kelvin,
            thermostat_time_scale_fs=thermostat_time_scale_fs,
            seed=seed,
        )
        saves["n"] += 1
        if saves["n"] <= kills:
            raise RuntimeError("simulated timeout")

    monkeypatch.setattr(md_mod, "_save_checkpoint", killing_save)
    return saves


def test_run_nvt_md_resume_reproduces(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A rollout interrupted at several checkpoints and resumed reproduces an
    uninterrupted run frame-for-frame (atol 1e-14), with no duplicate or missing
    frames.
    """
    from matbench_discovery.md import run_nvt_md_resumable

    atoms = bulk("Cu", cubic=True) * (2, 2, 2)
    run_kwargs = {
        "temperature_kelvin": 300,
        "n_steps": 60,  # 60 // 10 + 1 = 7 frames (odd -> last frame skips checkpoint)
        "record_interval": 10,
        "checkpoint_every_n_frames": 2,
        "progress_interval": 0,
    }
    reference = run_nvt_md_resumable(
        atoms, EMT(), out_path=str(tmp_path / "ref.extxyz"), **run_kwargs
    )
    assert len(reference) == 7

    # simulate timeouts: raise after each of the first 3 checkpoint writes, then let
    # the run finish; the outer loop reruns (resuming) until it completes
    saves = _fail_save_checkpoint_after(monkeypatch)

    out_path = str(tmp_path / "out.extxyz")
    resumed = None
    for _attempt in range(8):  # bounded resume loop
        try:
            resumed = run_nvt_md_resumable(
                atoms, EMT(), out_path=out_path, **run_kwargs
            )
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
    from matbench_discovery.md import run_nvt_md_resumable

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
    from matbench_discovery import md_models

    assert isinstance(load_calculator("emt"), EMT)
    assert isinstance(load_calculator("emt", dtype="float32"), EMT)  # dtype ignored

    seen_dtype = ""

    def dtype_aware(device: str, dtype: str = "float64") -> EMT:  # noqa: ARG001
        nonlocal seen_dtype
        seen_dtype = dtype
        return EMT()

    monkeypatch.setitem(
        md_models.MD_MODELS, "mace_mp_0", md_models.MdModel(dtype_aware)
    )
    load_calculator("mace_mp_0", device="cpu", dtype="float32")
    assert seen_dtype == "float32"

    with pytest.raises(ValueError, match="Unknown model_key"):
        load_calculator("does-not-exist")


def test_run_md_cli_write_yaml_skips_non_submission_model(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`--write-yaml` for a debug model with no Model enum entry (emt) must skip the
    YAML write gracefully instead of crashing on Model.from_ref.
    """
    import sys

    from models import run_md

    # avoid the real calculator + rollout pipeline; return metrics calc_md_metrics reads
    monkeypatch.setattr(run_md, "load_calculator", lambda *_a, **_k: EMT())
    monkeypatch.setattr(
        run_md,
        "run_md_benchmark",
        lambda **_kwargs: pd.DataFrame({"rdf_error": [1.0], "vdos_error": [2.0]}),
    )
    monkeypatch.setattr(sys, "argv", ["run_md", "--model", "emt", "--write-yaml"])

    assert run_md.main() == 0  # would raise ValueError without the non-enum guard


def test_run_md_cli_rejects_partial_write_yaml(monkeypatch: pytest.MonkeyPatch) -> None:
    """--write-yaml on a --systems subset must error before the rollout pipeline."""
    import sys

    from models import run_md

    argv = ["run_md", "--model", "mace_mp_0", "--write-yaml", "--systems", "bulkAu"]
    monkeypatch.setattr(sys, "argv", argv)
    with pytest.raises(SystemExit):  # parser.error exits before the rollout pipeline
        run_md.main()


def test_md_evals_skips_incomplete_coverage(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """scripts/evals/md.py must not write model YAML from incomplete coverage (only one
    of the three reference systems has a per-system CSV).
    """
    from scripts.evals import md as eval_md

    model = Model.mace_mp_0
    monkeypatch.setattr(eval_md, "ROOT", str(tmp_path))
    monkeypatch.setattr(eval_md.cli_args, "models", [model])
    ref_dir = tmp_path / "ref"
    for name in ("sysA", "sysB", "sysC"):
        (ref_dir / name).mkdir(parents=True)
    monkeypatch.setattr(
        eval_md, "default_md_reference_paths", lambda: (str(ref_dir), "settings.csv")
    )

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


def test_md_model_uv_run_cmd() -> None:
    """uv_run_cmd interleaves --with/--find-links around the script and its args."""
    cmd = MD_MODELS["orb_v3"].uv_run_cmd("models/run_md.py", "--model", "orb_v3")
    assert cmd[:2] == ["uv", "run"]
    assert cmd[-3:] == ["models/run_md.py", "--model", "orb_v3"]
    deps = [cmd[idx + 1] for idx, tok in enumerate(cmd) if tok == "--with"]
    assert deps == list(MD_MODELS["orb_v3"].deps)
    # a model with find_links must emit a --find-links token per URL
    fairchem = MD_MODELS["eqv2_s_dens_mp"]
    links_cmd = fairchem.uv_run_cmd("models/run_md.py", "--model", "eqv2_s_dens_mp")
    links = [
        links_cmd[idx + 1] for idx, tok in enumerate(links_cmd) if tok == "--find-links"
    ]
    assert links == list(fairchem.find_links)


@pytest.mark.parametrize("model_key", [key for key in MD_MODELS if key != "emt"])
def test_md_model_keys_are_valid_enum_members(model_key: str) -> None:
    """Every registered model (besides the emt debug entry) must map to a Model so
    metrics can be written to the right YAML.
    """
    assert Model.from_ref(model_key).name == model_key
