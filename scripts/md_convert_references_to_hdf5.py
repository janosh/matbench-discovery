"""Convert CFPMD-26 reference trajectories from extxyz(.xz) to HDF5 (one ``traj.h5``
per system) so the MD benchmark loads them in milliseconds instead of the minutes-to-
hours ASE's pure-Python extxyz parser needs for tens of thousands of frames.

Run once per dataset (the slow parse happens here, not in every benchmark job):

    python scripts/md_convert_references_to_hdf5.py            # all systems
    python scripts/md_convert_references_to_hdf5.py --systems bulkAu_1500K_Kapil
    python scripts/md_convert_references_to_hdf5.py --overwrite
"""

import argparse
import os
import sys
import time

from matbench_discovery.md import (
    default_md_reference_paths,
    find_reference_trajectory,
    read_trajectory,
)
from matbench_discovery.trajectory import Trajectory


def main() -> int:
    """Convert each reference system's extxyz trajectory to a sibling traj.h5."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-dir", help="Defaults to the CFPMD-26 reference dir")
    parser.add_argument("--systems", nargs="*", help="Subset of system dir names")
    parser.add_argument(
        "--overwrite", action="store_true", help="Re-convert even if traj.h5 exists"
    )
    args = parser.parse_args()

    ref_dir = args.ref_dir or default_md_reference_paths()[0]
    system_dirs = sorted(
        entry.name
        for entry in os.scandir(ref_dir)
        if entry.is_dir() and (not args.systems or entry.name in args.systems)
    )
    if not system_dirs:
        parser.error(f"No system directories found under {ref_dir!r}")

    n_converted = 0
    for system_name in system_dirs:
        h5_path = f"{ref_dir}/{system_name}/traj.h5"
        if os.path.isfile(h5_path) and not args.overwrite:
            print(f"skip {system_name}: {h5_path} exists (use --overwrite)")
            continue
        src = find_reference_trajectory(ref_dir, system_name)
        start = time.perf_counter()
        trajectory = Trajectory.from_ase(read_trajectory(src))
        trajectory.write_hdf5(h5_path)
        src_mb = os.path.getsize(src) / 1e6
        h5_mb = os.path.getsize(h5_path) / 1e6
        print(
            f"{system_name}: {trajectory.n_frames} frames, {os.path.basename(src)} "
            f"{src_mb:.1f} MB -> traj.h5 {h5_mb:.1f} MB in "
            f"{time.perf_counter() - start:.1f}s"
        )
        n_converted += 1

    print(f"\nConverted {n_converted}/{len(system_dirs)} systems under {ref_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
