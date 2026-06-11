"""Tests for shared molecular dynamics helpers."""

from pathlib import Path

import ase.io
from ase import Atoms

from matbench_discovery.md import (
    IpiRunConfig,
    MdRunConfig,
    extract_temperature_from_name,
    find_trajectory_files,
    infer_temperature,
    read_trajectory,
    should_skip_trajectory,
    write_ipi_molecular_crystal_workflows,
)


def test_extract_temperature_from_name() -> None:
    """Temperatures should be parsed from common system directory names."""
    assert extract_temperature_from_name("MAPbBr3_300K_Ivor_VASP") == 300
    assert extract_temperature_from_name("sample_295.5K") == 295.5
    assert extract_temperature_from_name("no-temperature") is None


def test_infer_temperature_prefers_atoms_info(tmp_path: Path) -> None:
    """Atoms metadata should provide temperature when present."""
    atoms = Atoms("Si", positions=[(0, 0, 0)], cell=[5, 5, 5], pbc=True)
    atoms.info["temperature_K"] = "350"
    assert infer_temperature(tmp_path / "sample_300K" / "traj.extxyz", atoms) == 350


def test_find_read_and_skip_trajectory_files(tmp_path: Path) -> None:
    """Trajectory helpers should find, read and skip expected files."""
    system_dir = tmp_path / "bulkSi_300K"
    system_dir.mkdir()
    traj_path = system_dir / "traj.extxyz"
    ase.io.write(
        traj_path, Atoms("Si", positions=[(0, 0, 0)], cell=[5, 5, 5], pbc=True)
    )

    skipped_dir = tmp_path / "H_100K"
    skipped_dir.mkdir()
    skipped_path = skipped_dir / "traj.extxyz"
    ase.io.write(
        skipped_path, Atoms("H", positions=[(0, 0, 0)], cell=[5, 5, 5], pbc=True)
    )

    files = find_trajectory_files(tmp_path)
    assert files == [skipped_path, traj_path]
    assert len(read_trajectory(traj_path)) == 1
    assert infer_temperature(traj_path) == 300

    config = MdRunConfig()
    assert should_skip_trajectory(skipped_path, config)
    assert not should_skip_trajectory(traj_path, config)


def test_write_ipi_molecular_crystal_workflows(tmp_path: Path) -> None:
    """i-PI workflow generation should write the four expected run files."""
    input_dir = tmp_path / "inputs"
    system_dir = input_dir / "anthracene_293K_Sharma_S"
    system_dir.mkdir(parents=True)
    traj_path = system_dir / "traj.extxyz"
    ase.io.write(
        traj_path,
        Atoms(
            "C",
            positions=[(0, 0, 0)],
            cell=[5, 5, 5],
            pbc=True,
        ),
    )
    ase.io.write(
        system_dir / "anthracene_trajectory.extxyz",
        Atoms("H", positions=[(0, 0, 0)], cell=[5, 5, 5], pbc=True),
    )

    output_dir = tmp_path / "ipi"
    df_results = write_ipi_molecular_crystal_workflows(
        model_key="mace-mp-0",
        input_dir=input_dir,
        output_dir=output_dir,
        config=IpiRunConfig(
            total_steps=10,
            total_time=10,
            timestep_fs=0.5,
            thermostat_tau_fs=10,
            thermostat_tau_fs_by_system={"anthracene": 500},
            cuda_visible_devices="1",
            overwrite=True,
        ),
        system_names=["anthracene_293K_Sharma_S"],
    )

    assert df_results["status"].tolist() == ["completed"]
    assert df_results["input_file"].tolist() == [str(traj_path)]
    run_dir = output_dir / "anthracene_293K_Sharma_S" / "mace-mp-0"
    assert len(read_trajectory(run_dir / "init.xyz")) == 1
    assert (run_dir / "input.xml").is_file()
    assert (run_dir / "run-ase.py").is_file()
    assert (run_dir / "submit.sh").is_file()

    input_xml = (run_dir / "input.xml").read_text(encoding="utf-8")
    assert "<total_steps>10</total_steps>" in input_xml
    assert "<temperature units='kelvin'> 293 </temperature>" in input_xml
    assert "500" in input_xml
    assert "driver_anthracene" in input_xml

    run_ase = (run_dir / "run-ase.py").read_text(encoding="utf-8")
    assert "mace_mp(" in run_ase
    assert 'os.getenv("MACE_MODEL", "medium")' in run_ase
    assert 'os.getenv("IPI_SOCKET", "driver_anthracene")' in run_ase

    submit = (run_dir / "submit.sh").read_text(encoding="utf-8")
    assert "#SBATCH" not in submit
    assert "source " not in submit
    assert "activation" not in submit
    assert "export CUDA_VISIBLE_DEVICES=1" in submit
    assert "export IPI_SOCKET=driver_anthracene" in submit
    assert "conda activate" not in submit
    assert "uv run i-pi input.xml" in submit
    assert "uv run python3 run-ase.py" in submit

    submit_all = (output_dir / "submit_all.sh").read_text(encoding="utf-8")
    assert "nohup bash submit.sh" in submit_all
    assert "sequential" in submit_all
