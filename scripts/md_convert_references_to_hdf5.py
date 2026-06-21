"""Pack the CFPMD-26 reference trajectories (per-system extxyz.xz + a settings CSV)
into a single multi-group HDF5: one group per system holding the trajectory arrays plus
``dt_fs`` (saved-frame interval) and ``temperature_kelvin`` as group attributes. This is
a one-time maintainer step run before hosting the reference dataset; the benchmark then
reads the HDF5 directly (no per-job extxyz parsing, no settings CSV).

    python scripts/md_convert_references_to_hdf5.py \
        --ref-dir ~/data/cfpmd-26/reference_AIMD_trajectories \
        --settings-csv ~/data/cfpmd-26/reference_AIMD_timestep_and_stride.csv \
        --out data/md/2026-06-12-cfpmd-26-aimd-reference-md-trajectories.h5

Parsing uses ``Trajectory.from_extxyz`` (bulk numpy/pandas reader), which converts even
the largest 20k-70k-frame references in seconds-to-minutes rather than the hours ASE's
per-frame extxyz parser needs, so a single process handles all systems sequentially.
Only the resulting HDF5 needs to be hosted/distributed.
"""

import argparse
import os
import time
from glob import glob
from typing import TYPE_CHECKING

import pandas as pd

from matbench_discovery.md import write_reference_h5
from matbench_discovery.trajectory import Trajectory

if TYPE_CHECKING:
    from collections.abc import Iterator


def load_settings(csv_path: str) -> dict[str, tuple[float, float]]:
    """Map '<System>_<temp>K' keys to (dt_fs, temperature_kelvin) from a CFPMD settings
    CSV with System, temperature, stride and dt columns.

    dt is the AIMD *integration* timestep (fs) and stride the number of MD steps
    between saved frames, so the interval between reference frames is ``dt * stride``
    (e.g. dt=1, stride=5 -> 5 fs). dt_fs must be this saved-frame interval: it sets the
    VDOS frequency axis and the reference/prediction time-matching, so dropping the
    stride factor (wrong by up to 5x here) corrupts the VDOS and pressure metrics.
    """
    return {
        f"{row['System']}_{float(row['temperature']):g}K": (
            float(row["dt"]) * float(row["stride"]),
            float(row["temperature"]),
        )
        for _, row in pd.read_csv(csv_path).iterrows()
    }


def resolve_settings(
    system_name: str, settings: dict[str, tuple[float, float]]
) -> tuple[float, float] | None:
    """(dt_fs, temperature_kelvin) for a trajectory dir name. Dir names append suffixes
    like '_Kapil' or '-Artrith_VASP' to the '<System>_<temp>K' settings key (delimiter
    '_' for most systems, '-' for a few); pick the longest delimiter-aware prefix match
    so a name containing another '<...>K' token doesn't resolve to a shorter key.
    """
    matches = {
        key: value
        for key, value in settings.items()
        if system_name == key
        or (system_name.startswith(key) and system_name[len(key)] in "_-")
    }
    return matches[max(matches, key=len)] if matches else None


def find_reference_trajectory(ref_dir: str, system_name: str) -> str:
    """Path to the single traj.extxyz (or ASE-transparent .gz/.xz/.bz2) under
    ref_dir/system_name.
    """
    traj_files = sorted(glob(f"{ref_dir}/{system_name}/traj.*xyz*"))
    if len(traj_files) != 1:
        raise ValueError(
            f"Expected 1 reference trajectory in {ref_dir}/{system_name}, "
            f"got {traj_files}"
        )
    return traj_files[0]


def main() -> int:
    """Pack per-system extxyz references + settings CSV into one multi-group HDF5."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-dir", required=True, help="Raw per-system extxyz dir")
    parser.add_argument(
        "--settings-csv", required=True, help="Per-system dt + temperature CSV"
    )
    parser.add_argument("--out", required=True, help="Output reference HDF5 path")
    parser.add_argument("--systems", nargs="*", help="Subset of system dir names")
    args = parser.parse_args()

    settings = load_settings(args.settings_csv)
    system_dirs = sorted(
        entry.name
        for entry in os.scandir(args.ref_dir)
        if entry.is_dir() and (not args.systems or entry.name in args.systems)
    )
    if not system_dirs:
        parser.error(f"No system directories found under {args.ref_dir!r}")

    # resolve settings for every system up front (cheap) so a missing CSV row fails
    # before the expensive parse and before the output file is opened
    resolved_raw = {name: resolve_settings(name, settings) for name in system_dirs}
    if missing := [name for name, value in resolved_raw.items() if value is None]:
        parser.error(f"No settings row in {args.settings_csv} for {missing}")
    resolved = {name: value for name, value in resolved_raw.items() if value}

    def entries() -> "Iterator[tuple[str, Trajectory, float, float]]":
        """Parse one system at a time so write_reference_h5 streams it to disk and
        only the single largest trajectory is ever held in memory.
        """
        for system_name in system_dirs:
            dt_fs, temperature_kelvin = resolved[system_name]
            src = find_reference_trajectory(args.ref_dir, system_name)
            start = time.perf_counter()
            trajectory = Trajectory.from_extxyz(src)
            print(
                f"{system_name}: {trajectory.n_frames} frames, dt={dt_fs} fs, "
                f"T={temperature_kelvin} K, parsed {os.path.basename(src)} "
                f"({os.path.getsize(src) / 1e6:.1f} MB) in "
                f"{time.perf_counter() - start:.1f}s",
                flush=True,
            )
            yield system_name, trajectory, dt_fs, temperature_kelvin

    write_reference_h5(args.out, entries())
    out_mb = os.path.getsize(args.out) / 1e6
    print(f"\nWrote {len(system_dirs)} systems to {args.out} ({out_mb:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
