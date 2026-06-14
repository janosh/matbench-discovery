"""Tests for molecular dynamics rollout helpers."""

from pathlib import Path

import ase.io
import numpy as np
import pandas as pd
import pytest
from ase.build import bulk
from ase.calculators.emt import EMT

from matbench_discovery.enums import Model
from matbench_discovery.md import (
    extract_temperature,
    load_frame_intervals,
    read_trajectory,
    resolve_frame_interval,
    run_md_benchmark,
    run_nvt_md,
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
        "bulkAg_600K_extra,600,1,9,32\n"
        "bulkCuAu,500,5,2,32\n"
    )
    intervals = load_frame_intervals(str(csv_path))
    assert intervals == {
        "bulkAg_600K": 1.0,
        "anthracene_293K": 0.5,
        "bulkAg_600K_extra_600K": 9.0,
        "bulkCuAu_500K": 2.0,
    }
    assert resolve_frame_interval("bulkAg_600K_Kapil", intervals) == 1.0
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


def test_run_md_benchmark_no_matching_systems(tmp_path: Path) -> None:
    """A system filter matching no directory raises a clear error rather than a
    confusing empty-DataFrame KeyError downstream.
    """
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


def test_find_per_system_csvs_prefix_safe(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Per-system CSV discovery should avoid model-prefix and combined-file matches."""
    from scripts.evals import md as eval_md

    monkeypatch.setattr(eval_md, "ROOT", str(tmp_path))
    md_dir_old = tmp_path / "models" / "orb" / "2026-06-13-md-nvt"
    md_dir_new = tmp_path / "models" / "orb" / "2026-06-14-md-nvt"
    for directory in (md_dir_old, md_dir_new):
        directory.mkdir(parents=True)

    expected_paths = [
        f"{md_dir_old}/orb_v2-md-metrics-bulkAg_600K.csv.gz",
        f"{md_dir_new}/orb_v2-md-metrics-bulkCu_600K.csv.gz",
    ]
    ignored_paths = [
        f"{md_dir_new}/orb_v2-md-metrics.csv.gz",  # combined file, no system suffix
        f"{md_dir_new}/orb_v2_mptrj-md-metrics-bulkAg_600K.csv.gz",  # prefix collision
    ]
    for csv_path in (*expected_paths, *ignored_paths):
        Path(csv_path).write_text("system,rdf_error\nx,1\n")

    assert eval_md.find_per_system_csvs(Model.orb_v2) == expected_paths


def test_load_calculator() -> None:
    """Emt loads with no extra deps; unknown keys raise a helpful error."""
    assert isinstance(load_calculator("emt"), EMT)
    with pytest.raises(ValueError, match="Unknown model_key"):
        load_calculator("does-not-exist")


def test_run_md_cli_write_yaml_skips_non_submission_model(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`--write-yaml` for a debug model with no Model enum entry (emt) must skip the
    YAML write gracefully instead of crashing on Model.from_ref.
    """
    import importlib.util
    import sys

    from matbench_discovery import ROOT

    spec = importlib.util.spec_from_file_location("run_md", f"{ROOT}/models/run_md.py")
    assert spec is not None
    assert spec.loader is not None
    run_md = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(run_md)

    # avoid the real calculator + rollout pipeline; return metrics calc_md_metrics reads
    monkeypatch.setattr(run_md, "load_calculator", lambda *_a, **_k: EMT())
    monkeypatch.setattr(
        run_md,
        "run_md_benchmark",
        lambda **_kwargs: pd.DataFrame({"rdf_error": [1.0], "vdos_error": [2.0]}),
    )
    monkeypatch.setattr(sys, "argv", ["run_md", "--model", "emt", "--write-yaml"])

    assert run_md.main() == 0  # would raise ValueError without the non-enum guard


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
