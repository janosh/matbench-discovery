"""Dense-array MD trajectory container with HDF5 I/O.

Decouples the MD benchmark from ASE/extxyz: reference trajectories of tens of
thousands of frames load from a chunked, compressed HDF5 file orders of magnitude
faster than the minutes-to-hours ASE's pure-Python extxyz parser needs to build one
``Atoms`` object per frame (measured ~4000x on a 625-frame reference). A full read
still returns the arrays in memory; strided/partial reads only touch the requested
chunks. The metric functions consume this columnar representation directly; ASE is
only used for the MD rollout (integrator + calculators) and interop (``to_ase`` /
``from_ase``).
"""

import os
from dataclasses import dataclass
from typing import TYPE_CHECKING

import h5py
import numpy as np

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

    from ase import Atoms

# bump when the on-disk layout changes incompatibly; readers fail closed on mismatch
TRAJECTORY_SCHEMA = 1
# per-frame array fields (None when absent); atomic_numbers/pbc are frame-independent
_FRAME_FIELDS = ("positions", "cell", "energy", "forces", "stress", "md_step")


@dataclass(frozen=True)
class Trajectory:
    """An MD trajectory as stacked numpy arrays (frame axis first).

    Args:
        atomic_numbers: (n_atoms,) int atomic numbers, constant across frames.
        positions: (n_frames, n_atoms, 3) Cartesian positions in Angstrom.
        cell: (n_frames, 3, 3) row-vector cells in Angstrom (per frame, so NPT works).
        pbc: (3,) bool periodic-boundary flags, constant across frames.
        energy: (n_frames,) potential energies in eV, or None.
        forces: (n_frames, n_atoms, 3) forces in eV/Angstrom, or None.
        stress: (n_frames, 6) Voigt stress [xx, yy, zz, yz, xz, xy] in eV/Å³, or None.
        md_step: (n_frames,) integer MD step counters, or None.
    """

    atomic_numbers: np.ndarray
    positions: np.ndarray
    cell: np.ndarray
    pbc: np.ndarray
    energy: "np.ndarray | None" = None
    forces: "np.ndarray | None" = None
    stress: "np.ndarray | None" = None
    md_step: "np.ndarray | None" = None

    def __post_init__(self) -> None:
        """Validate array shapes are mutually consistent."""
        n_frames, n_atoms = self.positions.shape[0], len(self.atomic_numbers)
        expected = {
            "positions": (n_frames, n_atoms, 3),
            "cell": (n_frames, 3, 3),
            "forces": (n_frames, n_atoms, 3),
            "energy": (n_frames,),
            "stress": (n_frames, 6),
            "md_step": (n_frames,),
        }
        for name, shape in expected.items():
            arr = getattr(self, name)
            if arr is not None and arr.shape != shape:
                raise ValueError(f"{name} has shape {arr.shape}, expected {shape}")
        if self.pbc.shape != (3,):
            raise ValueError(f"pbc has shape {self.pbc.shape}, expected (3,)")

    def __len__(self) -> int:
        """Number of frames."""
        return self.n_frames

    @property
    def n_frames(self) -> int:
        """Number of recorded frames."""
        return self.positions.shape[0]

    @property
    def n_atoms(self) -> int:
        """Number of atoms per frame."""
        return len(self.atomic_numbers)

    def __getitem__(self, frames: "slice | int") -> "Trajectory":
        """Sub-trajectory over a frame slice (int selects a single-frame trajectory,
        supporting Python's negative indexing).
        """
        if isinstance(frames, int):
            idx = frames + self.n_frames if frames < 0 else frames
            if not 0 <= idx < self.n_frames:
                raise IndexError(f"frame {frames} out of range for {self.n_frames}")
            frames = slice(idx, idx + 1)
        sliced = {
            name: arr[frames]
            for name in _FRAME_FIELDS
            if (arr := getattr(self, name)) is not None
        }
        return Trajectory(atomic_numbers=self.atomic_numbers, pbc=self.pbc, **sliced)

    def frame_as_atoms(self, idx: int) -> "Atoms":
        """The idx-th frame as an ASE ``Atoms`` (with results attached if present)."""
        return self[idx].to_ase()[0]

    def to_ase(self) -> "list[Atoms]":
        """Convert all frames to ASE ``Atoms`` with a SinglePointCalculator carrying
        energy/forces/stress when available. Used for interop and to feed calculators
        (energy/force RMSE); not on the hot read path.
        """
        from ase import Atoms
        from ase.calculators.singlepoint import SinglePointCalculator

        out: list[Atoms] = []
        for idx in range(self.n_frames):
            atoms = Atoms(
                numbers=self.atomic_numbers,
                positions=self.positions[idx],
                cell=self.cell[idx],
                pbc=self.pbc,
            )
            results = {
                key: arr[idx]
                for key in ("energy", "forces", "stress")
                if (arr := getattr(self, key)) is not None
            }
            if results:
                atoms.calc = SinglePointCalculator(atoms, **results)
            if self.md_step is not None:
                atoms.info["md_step"] = int(self.md_step[idx])
            out.append(atoms)
        return out

    @classmethod
    def from_ase(cls, frames: "Sequence[Atoms]") -> "Trajectory":
        """Build a Trajectory from ASE frames, extracting energy/forces/stress/md_step
        where present (e.g. from a SinglePointCalculator on extxyz frames).
        """
        from ase.calculators.calculator import PropertyNotImplementedError

        if len(frames) == 0:
            raise ValueError("Cannot build a Trajectory from zero frames")
        first = frames[0]
        n_atoms = len(first)
        for idx, atoms in enumerate(frames):
            if len(atoms) != n_atoms:
                raise ValueError(
                    f"Inconsistent atom counts: frame {idx} has {len(atoms)}, "
                    f"expected {n_atoms}"
                )
            if not np.array_equal(atoms.numbers, first.numbers):
                raise ValueError(f"frame {idx} chemical species differ from frame 0")
            if not np.array_equal(atoms.pbc, first.pbc):
                raise ValueError(f"frame {idx} pbc differs from frame 0")
        positions = np.array([atoms.positions for atoms in frames])
        cell = np.array([atoms.cell.array for atoms in frames])

        def stack(getter: "Callable[[Atoms], np.ndarray]") -> "np.ndarray | None":
            """Stack a per-frame quantity. Returns None if no frame has it; raises on
            mixed availability (so a partially-labeled trajectory fails loud instead
            of silently dropping the field) and on genuine calculator backend errors.
            """
            values: list[np.ndarray | None] = []
            for atoms in frames:
                if atoms.calc is None:
                    values.append(None)
                    continue
                try:
                    values.append(getter(atoms))
                except PropertyNotImplementedError:
                    values.append(None)  # property simply not stored on this frame
            present = [value is not None for value in values]
            if not any(present):
                return None
            if not all(present):
                raise ValueError(
                    f"frame {present.index(False)} missing a property "
                    "that other frames provide"
                )
            return np.array(values)

        md_steps = [atoms.info.get("md_step") for atoms in frames]
        has_steps = all(step is not None for step in md_steps)
        return cls(
            atomic_numbers=first.numbers.copy(),
            positions=positions,
            cell=cell,
            pbc=np.asarray(first.pbc, dtype=bool),
            energy=stack(lambda atoms: atoms.get_potential_energy()),
            forces=stack(lambda atoms: atoms.get_forces()),
            stress=stack(lambda atoms: atoms.get_stress(voigt=True)),
            md_step=np.array(md_steps) if has_steps else None,
        )

    def write_hdf5(self, path: str) -> None:
        """Atomically write to a gzip-compressed, frame-chunked HDF5 file (tmp +
        os.replace). Frame-major chunking keeps strided/partial reads cheap; metadata
        (atomic_numbers, pbc, schema) lives alongside the per-frame datasets.
        """
        tmp_path = f"{path}.tmp"
        chunk_frames = min(self.n_frames, 512) or 1
        chunks_for = {
            "positions": (chunk_frames, self.n_atoms, 3),
            "cell": (chunk_frames, 3, 3),
            "forces": (chunk_frames, self.n_atoms, 3),
            "energy": (chunk_frames,),
            "stress": (chunk_frames, 6),
            "md_step": (chunk_frames,),
        }
        with h5py.File(tmp_path, "w") as file:
            file.attrs["schema"] = TRAJECTORY_SCHEMA
            file.create_dataset("atomic_numbers", data=self.atomic_numbers)
            file.create_dataset("pbc", data=self.pbc)
            for name in _FRAME_FIELDS:
                if (arr := getattr(self, name)) is not None:
                    file.create_dataset(
                        name, data=arr, chunks=chunks_for[name], compression="gzip"
                    )
        os.replace(tmp_path, path)

    @classmethod
    def read_hdf5(cls, path: str, *, frames: slice = slice(None)) -> "Trajectory":
        """Read (a slice of) frames from an HDF5 file written by ``write_hdf5``.

        Only the requested frames are read/decompressed, so strided or head reads of
        huge trajectories are cheap. Fails closed on a schema mismatch.
        """
        with h5py.File(path, "r") as file:
            schema = int(file.attrs.get("schema", -1))
            if schema != TRAJECTORY_SCHEMA:
                raise ValueError(
                    f"{path} has trajectory schema {schema}, expected "
                    f"{TRAJECTORY_SCHEMA}; regenerate it from the source trajectory"
                )
            kwargs = {
                name: file[name][frames] for name in _FRAME_FIELDS if name in file
            }
            return cls(
                atomic_numbers=file["atomic_numbers"][:],
                pbc=file["pbc"][:].astype(bool),
                **kwargs,
            )
