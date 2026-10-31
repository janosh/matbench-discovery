"""Tests for scripts/md_convert_references_to_hdf5.py."""

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


SYSTEM = "bulkCu_300K-Artrith_VASP"


def write_raw_reference(tmp_path: Path, dt_fs: float = 0.25) -> list[ase.Atoms]:
    """Write one short EMT Cu NVT trajectory + settings CSV under tmp_path."""
    frames = run_nvt_md(
        bulk("Cu", cubic=True) * (2, 2, 2),
        EMT(),
        temperature_kelvin=300,
        n_steps=8,
        record_interval=1,
    )
    sys_dir = tmp_path / "raw" / SYSTEM
    sys_dir.mkdir(parents=True)
    ase.io.write(sys_dir / "traj.extxyz", frames)
    settings_csv = tmp_path / "settings.csv"
    settings_csv.write_text(f"System,temperature,stride,dt\nbulkCu,300,4,{dt_fs}\n")
    return frames


def packer_argv(tmp_path: Path, *extra: str) -> list[str]:
    """CLI args packing tmp_path/raw with tmp_path/settings.csv into ref.h5."""
    return [
        *("--ref-dir", str(tmp_path / "raw")),
        *("--settings-csv", str(tmp_path / "settings.csv")),
        *("--out", str(tmp_path / "ref.h5")),
        *extra,
    ]


def test_md_convert_packs_reference(tmp_path: Path) -> None:
    """The packer turns raw per-system extxyz + settings CSV into reference HDF5."""
    frames = write_raw_reference(tmp_path)
    system = SYSTEM
    out_h5 = str(tmp_path / "ref.h5")
    assert PACKER.main(packer_argv(tmp_path)) == 0

    assert list_reference_systems(out_h5) == [system]
    traj, dt_fs, temperature = read_reference_trajectory(out_h5, system)
    assert dt_fs == 0.25  # dt = 0.25 fs (stride 4 not folded into saved-frame interval)
    assert temperature == 300.0
    assert traj.n_frames == len(frames)
    # extxyz round-trips positions to print precision
    np.testing.assert_allclose(
        traj.positions[0], frames[0].positions, rtol=0, atol=1e-6
    )


@pytest.mark.parametrize(
    ("dt_fs", "extra_args", "error", "match"),
    [
        # a typo'd --systems name must not silently pack fewer systems
        (0.25, ("--systems", SYSTEM, "bulkCu_typo"), SystemExit, "2"),
        # 10x too large dt_fs underestimates equipartition T by 100x
        (2.5, (), ValueError, "dt_fs=2.5 fs .* is likely wrong"),
    ],
)
def test_md_convert_rejects_bad_inputs(
    tmp_path: Path,
    dt_fs: float,
    extra_args: tuple[str, ...],
    error: type[Exception],
    match: str,
) -> None:
    """Unknown --systems and wrong dt_fs fail instead of writing the HDF5."""
    write_raw_reference(tmp_path, dt_fs=dt_fs)
    with pytest.raises(error, match=match):
        PACKER.main(packer_argv(tmp_path, *extra_args))
    assert not (tmp_path / "ref.h5").exists()
