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
    energy: "np.ndarray | None" = None
    forces: "np.ndarray | None" = None
    stress: "np.ndarray | None" = None
    md_step: "np.ndarray | None" = None

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
        present = [step is not None for step in md_steps]
        # match stack(): all present -> keep, none present -> drop, partial -> fail loud
        # rather than silently dropping the field
        if any(present) and not all(present):
            frame_idx = present.index(False)
            raise ValueError(f"{frame_idx=} missing md_step that other frames have")
        return cls(
            atomic_numbers=first.numbers.copy(),
            positions=positions,
            cell=cell,
            pbc=np.asarray(first.pbc, dtype=bool),
            energy=stack(lambda atoms: atoms.get_potential_energy()),
            forces=stack(lambda atoms: atoms.get_forces()),
            stress=stack(lambda atoms: atoms.get_stress(voigt=True)),
            md_step=np.array(md_steps) if all(present) else None,
        )

    @classmethod
    def from_extxyz(cls, path: str) -> "Trajectory":
        """Bulk-parse a (optionally xz/gz/bz2-compressed) extxyz trajectory into stacked
        arrays without building one ASE ``Atoms`` per frame.

        ASE's pure-Python extxyz reader rebuilds an ``Atoms`` (and re-parses the header)
        for every frame, which is minutes-to-hours for the 10k-70k-frame CFPMD-26
        references (~0.5 frames/s measured). This reads the file once, slices the
        constant-stride atom blocks with numpy, and parses the numeric columns with
        pandas' C engine, reproducing ASE's cell/stress/energy conventions exactly
        (verified against ``from_ase(read_trajectory(...))`` in tests). It assumes a
        constant atom count and column layout across frames (true for these NVT
        references); a varying atom count fails loud on the frame-stride check.
        """
        import bz2
        import gzip
        import lzma
        import re
        from io import StringIO

        import pandas as pd
        from ase.symbols import symbols2numbers

        openers = {".xz": lzma.open, ".gz": gzip.open, ".bz2": bz2.open}
        opener = openers.get(os.path.splitext(path)[1], open)
        with opener(path, mode="rt") as file:
            lines = file.read().split("\n")
        while lines and not lines[-1].strip():  # tolerate trailing blank lines
            lines.pop()
        if not lines:
            raise ValueError(f"empty extxyz file {path!r}")

        n_atoms = int(lines[0].split()[0])
        stride = n_atoms + 2  # count line + comment line + n_atoms atom lines
        n_frames, remainder = divmod(len(lines), stride)
        if remainder or n_frames == 0:
            raise ValueError(
                f"{path!r}: {len(lines)} lines is not a multiple of frame stride "
                f"{stride} (n_atoms={n_atoms}); from_extxyz needs a constant atom count"
            )
        grid = np.asarray(lines, dtype=object).reshape(n_frames, stride)
        comments = grid[:, 1]
        atom_lines = grid[:, 2:].reshape(-1)  # frame-major (n_frames * n_atoms,)

        # frame-invariant metadata is parsed from frame 0 only (count line, Properties
        # header, species, pbc); validate every frame matches so differing later frames
        # can't silently mix systems into the stacked arrays
        counts = grid[:, 0].astype(int)
        if (count_differs := counts != n_atoms).any():
            bad = int(np.argmax(count_differs))
            raise ValueError(
                f"{path!r}: frame {bad} count line says {counts[bad]} atoms != "
                f"{n_atoms} (frame 0)"
            )

        def uniform_comment(pattern: str, label: str, default: str) -> str:
            """Frame 0's `pattern` capture, raising if any frame differs."""
            values = [
                m.group(1) if (m := re.search(pattern, c)) else default
                for c in comments
            ]
            for idx, value in enumerate(values):
                if value != values[0]:
                    raise ValueError(
                        f"{path!r}: frame {idx} {label} differs from frame 0"
                    )
            return values[0]

        # column layout from the Properties header (name:type:count triples)
        properties = uniform_comment(r"Properties=(\S+)", "Properties header", "")
        if not properties:
            raise ValueError(f"{path!r}: first comment line has no Properties=...")
        tokens = properties.split(":")
        col_start: dict[str, int] = {}
        col = 0
        for name, count in zip(tokens[::3], tokens[2::3], strict=True):
            col_start[name] = col
            col += int(count)
        for required in ("species", "pos"):
            if required not in col_start:
                raise ValueError(f"{path!r}: Properties is missing {required!r}")

        # bulk-parse every atom block in one C-engine call (\s+ skips leading space)
        table = pd.read_csv(
            StringIO("\n".join(atom_lines.tolist())),
            sep=r"\s+",
            header=None,
            names=range(col),
        )

        def grab(start: int) -> np.ndarray:
            """The 3 columns at ``start`` as (n_frames, n_atoms, 3) float."""
            block = table.iloc[:, start : start + 3].to_numpy(dtype=float)
            return block.reshape(n_frames, n_atoms, 3)

        positions = grab(col_start["pos"])
        forces = grab(col_start["forces"]) if "forces" in col_start else None
        species_col = table.iloc[:, col_start["species"]].to_numpy(dtype=str)
        species_all = species_col.reshape(n_frames, -1)
        if (species_differs := (species_all != species_all[0]).any(axis=1)).any():
            bad = int(np.argmax(species_differs))
            raise ValueError(f"{path!r}: frame {bad} species differ from frame 0")
        species = species_all[0]
        atomic_numbers = np.asarray(symbols2numbers(list(species)))

        # per-frame floats parsed from the comment lines
        def comment_field(
            pattern: str, width: int, *, drop_if_partial: bool = False
        ) -> "np.ndarray | None":
            """Flat (n_frames * width,) floats from a per-frame comment field, or None
            if no frame carries it. Mixed availability raises (matching from_ase) unless
            ``drop_if_partial``, which returns None instead (md_step's looser rule).
            """
            rgx = re.compile(pattern)
            found = [m.group(1) if (m := rgx.search(c)) else None for c in comments]
            n_missing = found.count(None)
            if n_missing == n_frames or (n_missing and drop_if_partial):
                return None  # absent on every frame, or partial for an optional field
            if n_missing:
                raise ValueError(
                    f"{path!r}: frame {found.index(None)} is missing a comment "
                    "field that other frames provide"
                )
            flat = np.fromstring(" ".join(v for v in found if v is not None), sep=" ")
            if flat.size != n_frames * width:
                raise ValueError(
                    f"{path!r}: expected {n_frames * width} values for {pattern!r}, "
                    f"got {flat.size}"
                )
            return flat

        lattices = comment_field(r'Lattice="([^"]*)"', 9)
        if lattices is None:
            raise ValueError(f"{path!r}: every frame must carry a Lattice")
        cell = lattices.reshape(n_frames, 3, 3)  # 9 row-major floats == ASE cell

        # (?<!\w): only a bare energy= key, never *_energy= (free_/kinetic_/total_...)
        energy = comment_field(r"(?<!\w)energy=(\S+)", 1)  # width 1 -> (n_frames,)

        stress = comment_field(r'\bstress="([^"]*)"', 9)
        if stress is not None:
            # ASE Fortran-reshapes the 9-vec, then Voigt-orders [xx,yy,zz,yz,xz,xy]
            stress = stress.reshape(n_frames, 9)[:, [0, 4, 8, 7, 6, 3]]

        steps = comment_field(r"\bmd_step=(\S+)", 1, drop_if_partial=True)
        md_step = steps.astype(int) if steps is not None else None

        pbc_str = uniform_comment(r'pbc="([^"]*)"', "pbc", "T T T")
        pbc = np.array([tok.lower().startswith("t") for tok in pbc_str.split()])

        return cls(
            atomic_numbers=atomic_numbers,
            positions=positions,
            cell=cell,
            pbc=pbc,
            energy=energy,
            forces=forces,
            stress=stress,
            md_step=md_step,
        )

    def write_to_h5_group(self, group: "h5py.Group") -> None:
        """Write this trajectory's arrays + metadata into an open HDF5 group.

        An ``h5py.File`` is also a ``Group``, so this backs both single-trajectory
        files (``write_hdf5``) and per-system groups in a multi-trajectory reference
        file. Frame-major chunking keeps strided/partial reads cheap.
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
        cls, group: "h5py.Group", *, frames: slice = slice(None)
    ) -> "Trajectory":
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

    def write_hdf5(self, path: str) -> None:
        """Atomically write to a gzip-compressed, frame-chunked HDF5 file (tmp +
        os.replace). Frame-major chunking keeps strided/partial reads cheap; metadata
        (atomic_numbers, pbc, schema) lives alongside the per-frame datasets.
        """
        tmp_path = f"{path}.tmp"
        with h5py.File(tmp_path, "w") as file:
            self.write_to_h5_group(file)
        os.replace(tmp_path, path)

    @classmethod
    def read_hdf5(cls, path: str, *, frames: slice = slice(None)) -> "Trajectory":
        """Read (a slice of) frames from an HDF5 file written by ``write_hdf5``."""
        with h5py.File(path, "r") as file:
            return cls.read_from_h5_group(file, frames=frames)
