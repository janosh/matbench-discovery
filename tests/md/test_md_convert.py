"""Tests for scripts/md_convert_references_to_hdf5.py."""

import sys
from pathlib import Path

import ase.io
import numpy as np
import pytest
from ase.build import bulk
from ase.calculators.emt import EMT

from matbench_discovery.md import (
    list_reference_systems,
    read_reference_trajectory,
    run_nvt_md,
)
from tests.utils import import_repo_script

PACKER = import_repo_script(
    "md_convert_references_to_hdf5", "scripts/md_convert_references_to_hdf5.py"
)


def test_md_convert_resolve_settings() -> None:
    """resolve_settings does suffix-aware longest-prefix dir-name to key matching."""
    settings = {
        "bulkAg_600K": (1.0, 600.0),
        "bulkAg_600K_extra_600K": (9.0, 600.0),
        "bulkCuAu_500K": (2.0, 500.0),
    }
    assert PACKER.resolve_settings("bulkAg_600K_Kapil", settings) == (1.0, 600.0)
    assert PACKER.resolve_settings("bulkAg_600K_extra_600K_Kapil", settings) == (
        9.0,
        600.0,
    )
    assert PACKER.resolve_settings("bulkAg_600K2_Kapil", settings) is None
    # '-' delimiter after the key must still resolve (caught a prod bug)
    assert PACKER.resolve_settings("bulkCuAu_500K-Artrith_VASP", settings) == (
        2.0,
        500.0,
    )


def test_md_convert_load_settings(tmp_path: Path) -> None:
    """load_settings reads dt_fs = the saved-frame interval = the CSV's dt column.
    The reference extxyz stores every integration step, so the stride column is NOT
    folded in (verified by equipartition; dt * stride underestimates temperature by
    stride**2 and collapses the vDOS frequency axis).
    """
    csv = tmp_path / "settings.csv"
    csv.write_text(
        "System,temperature,stride,dt\nTiSe2,400,5,1\nanthracene,293,1,0.5\n"
    )
    settings = PACKER.load_settings(str(csv))
    assert settings["TiSe2_400K"] == (1.0, 400.0)  # dt = 1 fs (stride 5 not folded in)
    assert settings["anthracene_293K"] == (0.5, 293.0)  # dt = 0.5 fs


def test_md_convert_packs_reference(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The packer turns raw per-system extxyz + settings CSV into reference HDF5."""
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
    assert PACKER.main() == 0

    assert list_reference_systems(out_h5) == [system]
    traj, dt_fs, temperature = read_reference_trajectory(out_h5, system)
    assert dt_fs == 0.25  # dt = 0.25 fs (stride 4 not folded into saved-frame interval)
    assert temperature == 300.0
    assert traj.n_frames == len(frames)
    # extxyz round-trips positions to print precision
    np.testing.assert_allclose(
        traj.positions[0], frames[0].positions, rtol=0, atol=1e-6
    )
