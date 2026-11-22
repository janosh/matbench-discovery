"""Dense-array MD trajectory container with HDF5 I/O.

Reference trajectories of tens of thousands of frames load from a chunked, compressed
HDF5 file orders of magnitude faster than parsing extxyz into one object per frame
(measured ~4000x on a 625-frame reference with ASE's reader). A full read still returns
the arrays in memory; strided/partial reads only touch the requested chunks. The metric
functions consume this columnar representation directly; per-frame
``kemio.io.Frame`` objects are only built for MD rollouts and interop
(``to_frames`` / ``from_frames``).
"""

from dataclasses import dataclass
from typing import TYPE_CHECKING

import h5py
import numpy as np

if TYPE_CHECKING:
    from collections.abc import Sequence

    from kemio import Structure
    from kemio.io import Frame

# bump when the on-disk layout changes incompatibly; readers fail closed on mismatch
TRAJECTORY_SCHEMA = 1
# per-frame array fields (None when absent); atomic_numbers/pbc are frame-independent
_FRAME_FIELDS = ("positions", "cell", "energy", "forces", "stress", "md_step")


# expected shape of each per-frame field for a leading-axis length (n_frames when
# validating, the chunk size when writing HDF5) and atom count
_frame_field_shapes = lambda n_lead, n_atoms: {  # noqa: E731
    "positions": (n_lead, n_atoms, 3),
    "cell": (n_lead, 3, 3),
    "forces": (n_lead, n_atoms, 3),
    "energy": (n_lead,),
    "stress": (n_lead, 6),
    "md_step": (n_lead,),
}


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
    energy: np.ndarray | None = None
    forces: np.ndarray | None = None
    stress: np.ndarray | None = None
    md_step: np.ndarray | None = None

    def __post_init__(self) -> None:
        """Validate array shapes are mutually consistent."""
        n_frames, n_atoms = self.positions.shape[0], len(self.atomic_numbers)
        expected = _frame_field_shapes(n_frames, n_atoms)
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

    def __getitem__(self, frames: slice | int) -> Trajectory:
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

    def frame_as_structure(self, idx: int) -> Structure:
        """The idx-th frame as a ``kemio.Structure`` (md_step in ``info``)."""
        return self[idx].to_frames()[0].structure

    def to_frames(self) -> list[Frame]:
        """Convert all frames to ``kemio.io.Frame`` objects carrying
        energy/forces/stress when available. Used for interop and to feed models
        (energy/force RMSE); not on the hot read path.
        """
        from kemio import Structure
        from kemio.elements import CHEMICAL_SYMBOLS
        from kemio.io import Frame
        from kemio.model import voigt_to_matrix

        symbols = [CHEMICAL_SYMBOLS[number] for number in self.atomic_numbers]
        out: list[Frame] = []
        for idx in range(self.n_frames):
            info = {} if self.md_step is None else {"md_step": int(self.md_step[idx])}
            structure = Structure(
                self.cell[idx],
                symbols,
                self.positions[idx],
                pbc=tuple(bool(flag) for flag in self.pbc),
                info=info,
            )
            out.append(
                Frame(
                    structure,
                    energy=None if self.energy is None else float(self.energy[idx]),
                    forces=None if self.forces is None else self.forces[idx],
                    stress=None
                    if self.stress is None
                    else voigt_to_matrix(self.stress[idx]),
                )
            )
        return out

    @classmethod
    def from_frames(cls, frames: Sequence[Frame]) -> Trajectory:
        """Build a Trajectory from frames (e.g. MD rollouts; files go through
        ``from_extxyz``) with ``kemio.interop.trajectory_arrays``, which raises
        ``NotStackableError`` (a ValueError) for no frames or frames that differ in
        atoms or pbc. Energy, forces, stress or md_step on only some frames raise
        ValueError.
        """
        from kemio.interop import trajectory_arrays
        from kemio.model import matrix_to_voigt

        arrays = trajectory_arrays(list(frames))
        md_step = arrays.info.get("md_step")
        fields = {
            "energy": arrays.energy,
            "forces": arrays.forces,
            "stress": arrays.stress,
            "md_step": md_step,
        }
        # trajectory_arrays fills a quantity missing in some frames with NaN
        if partial := [
            name
            for name, values in fields.items()
            if values is not None and np.isnan(values).any()
        ]:
            raise ValueError(f"only some frames carry {partial}")
        return cls(
            atomic_numbers=np.asarray(arrays.numbers, dtype=int),
            positions=arrays.positions,
            # the frames' lattices also for pbc=False (arrays.cells is None then)
            cell=np.array([frame.structure.lattice for frame in frames]),
            pbc=np.asarray(arrays.pbc, dtype=bool),
            energy=arrays.energy,
            forces=arrays.forces,
            stress=None if arrays.stress is None else matrix_to_voigt(arrays.stress),
            md_step=md_step,
        )

    @classmethod
    def from_extxyz(cls, path: str) -> Trajectory:
        """Read a (optionally .gz/.xz/.bz2) one-system extxyz trajectory with
        kemio's parallel Rust reader straight into stacked arrays (positions exactly
        as written). Frames that differ in atoms, pbc or columns, or that only partly
        carry energy, forces, stress or md_step, raise ValueError.
        """
        from kemio.io import read_extxyz_arrays
        from kemio.model import matrix_to_voigt

        arrays = read_extxyz_arrays(path)
        stress = arrays["stress"]
        fields = ("pbc", "positions", "cell", "energy", "forces")
        return cls(
            # kemio reads uint8 atomic numbers; keep the int64 of the other readers
            atomic_numbers=np.asarray(arrays["atomic_numbers"], dtype=int),
            **{name: arrays[name] for name in fields},
            md_step=arrays["md_step"],
            stress=None if stress is None else matrix_to_voigt(stress),
        )

    def write_to_h5_group(self, group: h5py.Group) -> None:
        """Write this trajectory's arrays + metadata into an open HDF5 group.

        An ``h5py.File`` is also a ``Group``, so this backs both single-trajectory
        files and per-system groups in a multi-trajectory reference file. Frame-major
        chunking keeps strided/partial reads cheap.
        """
        chunk_frames = min(self.n_frames, 512) or 1
        chunks_for = _frame_field_shapes(chunk_frames, self.n_atoms)
        group.attrs["schema"] = TRAJECTORY_SCHEMA
        group.create_dataset("atomic_numbers", data=self.atomic_numbers)
        group.create_dataset("pbc", data=self.pbc)
        for name in _FRAME_FIELDS:
            if (arr := getattr(self, name)) is not None:
                group.create_dataset(
                    name, data=arr, chunks=chunks_for[name], compression="gzip"
                )

    @classmethod
    def read_from_h5_group(
        cls, group: h5py.Group, *, frames: slice = slice(None)
    ) -> Trajectory:
        """Read (a slice of) frames from an HDF5 group written by ``write_to_h5_group``.

        Only the requested frames are read/decompressed, so strided or head reads of
        huge trajectories are cheap. Fails closed on a schema mismatch.
        """
        schema = int(group.attrs.get("schema", -1))
        if schema != TRAJECTORY_SCHEMA:
            raise ValueError(
                f"{group.name} has trajectory schema {schema}, expected "
                f"{TRAJECTORY_SCHEMA}; regenerate it from the source trajectory"
            )
        kwargs = {name: group[name][frames] for name in _FRAME_FIELDS if name in group}
        return cls(
            atomic_numbers=group["atomic_numbers"][:],
            pbc=group["pbc"][:].astype(bool),
            **kwargs,
        )
