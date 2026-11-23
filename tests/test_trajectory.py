"""Tests for the dense-array Trajectory container and its HDF5 / ASE interop."""

from pathlib import Path
from typing import Any

import ase.io
import h5py
import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from matbench_discovery.md import read_trajectory
from matbench_discovery.trajectory import TRAJECTORY_SCHEMA, Trajectory

rng = np.random.default_rng(seed=0)


def make_traj(
    n_frames: int = 6, n_atoms: int = 4, *, results: bool = True
) -> Trajectory:
    """Synthetic Trajectory with optional energy/forces/stress/md_step."""
    fields: dict[str, np.ndarray] = {
        "atomic_numbers": np.array([29] * n_atoms),
        "positions": rng.random((n_frames, n_atoms, 3)),
        "cell": np.broadcast_to(np.eye(3) * 5, (n_frames, 3, 3)).copy(),
        "pbc": np.array([True, True, True]),
    }
    if results:
        fields |= {
            "energy": rng.random(n_frames),
            "forces": rng.random((n_frames, n_atoms, 3)),
            "stress": rng.random((n_frames, 6)),
            "md_step": np.arange(n_frames) * 10,
        }
    return Trajectory(**fields)


def make_atoms(
    symbols: str = "Cu2",
    *,
    pbc: bool = True,
    energy: float | None = None,
    md_step: int | None = None,
) -> Atoms:
    """Minimal ASE Atoms (positions at the origin) in a 4 A cubic cell, optionally
    carrying an energy label or md_step.
    """
    atoms = Atoms(symbols, cell=np.eye(3) * 4, pbc=pbc)
    if energy is not None:
        atoms.calc = SinglePointCalculator(atoms, energy=energy)
    if md_step is not None:
        atoms.info["md_step"] = md_step
    return atoms


def extxyz_frame(
    *,
    count: int = 2,
    n_lines: int = 2,
    species: str = "Cu",
    props: str = "species:S:1:pos:R:3",
    pbc: str = "T T T",
    extra: str = "",
) -> str:
    """One extxyz frame (atom idx at (idx, idx, idx)) in a 4 A cubic cell."""
    header = f'Lattice="4 0 0 0 4 0 0 0 4" Properties={props} pbc="{pbc}"{extra}'
    atom_lines = "".join(f"{species} {idx} {idx} {idx}\n" for idx in range(n_lines))
    return f"{count}\n{header}\n{atom_lines}"


@pytest.mark.parametrize(
    ("field", "bad_shape"),
    [
        ("positions", (6, 3, 3)),  # n_atoms mismatch vs atomic_numbers
        ("cell", (6, 2, 3)),
        ("forces", (5, 4, 3)),  # n_frames mismatch
        ("energy", (5,)),
        ("stress", (6, 3)),
        ("md_step", (5,)),
    ],
)
def test_trajectory_validates_shapes(field: str, bad_shape: tuple[int, ...]) -> None:
    """__post_init__ rejects per-frame arrays with inconsistent shapes."""
    fields: dict[str, np.ndarray] = {
        "atomic_numbers": np.array([29] * 4),
        "positions": np.zeros((6, 4, 3)),
        "cell": np.zeros((6, 3, 3)),
        "pbc": np.array([True, True, True]),
        field: np.zeros(bad_shape),
    }
    with pytest.raises(ValueError, match="shape"):
        Trajectory(**fields)


def test_trajectory_len_props_and_slicing() -> None:
    """len/n_frames/n_atoms and indexing/slicing return consistent sub-trajectories."""
    traj = make_traj(n_frames=6, n_atoms=4)
    assert len(traj) == traj.n_frames == 6
    assert traj.n_atoms == 4

    sub = traj[1:4]
    assert sub.n_frames == 3
    assert traj.md_step is not None
    np.testing.assert_array_equal(sub.positions, traj.positions[1:4])
    np.testing.assert_array_equal(sub.md_step, traj.md_step[1:4])
    assert (sub.atomic_numbers == traj.atomic_numbers).all()

    single = traj[2]  # int selects a 1-frame trajectory
    assert single.n_frames == 1
    np.testing.assert_array_equal(single.positions[0], traj.positions[2])
    np.testing.assert_array_equal(traj[-1].positions[0], traj.positions[5])
    np.testing.assert_array_equal(traj.frame_as_atoms(-1).positions, traj.positions[5])
    with pytest.raises(IndexError):
        _ = traj[6]
    with pytest.raises(IndexError):
        _ = traj[-7]


@pytest.mark.parametrize("results", [True, False])
def test_trajectory_hdf5_roundtrip(tmp_path: Path, *, results: bool) -> None:
    """Full and strided HDF5 reads reproduce the written arrays exactly."""
    traj = make_traj(n_frames=20, n_atoms=5, results=results)
    path = str(tmp_path / "traj.h5")
    with h5py.File(path, "w") as file:
        traj.write_to_h5_group(file)

    with h5py.File(path, "r") as file:
        full = Trajectory.read_from_h5_group(file)
    np.testing.assert_array_equal(full.positions, traj.positions)
    np.testing.assert_array_equal(full.cell, traj.cell)
    assert (full.atomic_numbers == traj.atomic_numbers).all()
    assert (full.pbc == traj.pbc).all()
    for field in ("energy", "forces", "stress", "md_step"):
        if results:
            np.testing.assert_array_equal(getattr(full, field), getattr(traj, field))
        else:
            assert getattr(full, field) is None

    with h5py.File(path, "r+") as file:
        strided = Trajectory.read_from_h5_group(file, frames=slice(0, 20, 4))
        file.attrs["schema"] = TRAJECTORY_SCHEMA + 1
        with pytest.raises(ValueError, match="trajectory schema"):
            Trajectory.read_from_h5_group(file)
    assert strided.n_frames == 5
    np.testing.assert_array_equal(strided.positions, traj.positions[0:20:4])


def test_trajectory_ase_roundtrip() -> None:
    """from_ase extracts results; to_ase restores positions/cell/results/md_step."""
    traj = make_traj(n_frames=5, n_atoms=3)
    assert traj.forces is not None
    assert traj.energy is not None
    assert traj.stress is not None
    assert traj.md_step is not None
    atoms_list = traj.to_ase()
    assert len(atoms_list) == 5
    np.testing.assert_array_equal(atoms_list[0].get_forces(), traj.forces[0])
    assert atoms_list[2].info["md_step"] == traj.md_step[2]

    back = Trajectory.from_ase(atoms_list)
    assert back.energy is not None
    assert back.stress is not None
    np.testing.assert_allclose(back.positions, traj.positions, rtol=0, atol=1e-12)
    np.testing.assert_allclose(back.energy, traj.energy, rtol=0, atol=1e-12)
    np.testing.assert_allclose(back.stress, traj.stress, rtol=0, atol=1e-12)


def test_from_ase_missing_results() -> None:
    """ASE conversions omit calculators/results when frames carry no results."""
    traj = Trajectory.from_ase([make_atoms()])
    assert traj.energy is None
    assert traj.forces is None
    assert traj.stress is None
    assert traj.md_step is None

    atoms = make_traj(results=False).to_ase()[0]
    assert atoms.calc is None


@pytest.mark.parametrize(
    ("frames", "match"),
    [
        ([], "zero frames"),
        ([make_atoms("Cu2"), make_atoms("Cu3")], "Inconsistent atom counts"),
        ([make_atoms("Cu2"), make_atoms("He2")], "species differ"),
        ([make_atoms("Cu2"), make_atoms("Cu2", pbc=False)], "pbc differs"),
        # a property on some frames but not others fails loud (not silently dropped)
        ([make_atoms(energy=-1.0), make_atoms()], "frame 1 missing energy"),
        ([make_atoms(), make_atoms(md_step=0)], "frame 0 missing md_step"),
    ],
    ids=["empty", "ragged_count", "species", "pbc", "mixed_energy", "mixed_md_step"],
)
def test_from_ase_rejects_bad_input(frames: list[Atoms], match: str) -> None:
    """from_ase rejects empty input and frames inconsistent with frame 0."""
    with pytest.raises(ValueError, match=match):
        Trajectory.from_ase(frames)


@pytest.mark.parametrize("suffix", ["", ".xz", ".gz"])
def test_from_extxyz_matches_ase_parser(tmp_path: Path, suffix: str) -> None:
    """The fast bulk reader matches the ASE extxyz reader, compressed or not."""
    traj = make_traj(n_frames=5, n_atoms=4)
    assert traj.energy is not None
    # non-orthogonal cell so the parity check catches Lattice row/column-order or
    # transpose bugs that a diagonal cell would hide
    traj.cell[:] = [[5.0, 0.5, 0.3], [0.0, 5.2, 0.4], [0.0, 0.0, 4.7]]
    # out-of-cell positions catch wrapping (MD trajectories stay unwrapped)
    traj.positions[:, 0] += [12.0, -7.0, 30.0]
    path = str(tmp_path / f"traj.extxyz{suffix}")
    ase.io.write(path, traj.to_ase())

    fast = Trajectory.from_extxyz(path)
    slow = Trajectory.from_ase(read_trajectory(path))
    assert slow.energy is not None
    assert slow.n_frames == 5
    # extxyz is a text format, so positions/energy round-trip to its print precision
    np.testing.assert_allclose(slow.positions, traj.positions, rtol=0, atol=1e-6)
    np.testing.assert_allclose(slow.energy, traj.energy, rtol=0, atol=1e-6)
    np.testing.assert_array_equal(fast.atomic_numbers, slow.atomic_numbers)
    np.testing.assert_array_equal(fast.pbc, slow.pbc)
    for field in ("positions", "cell", "energy", "forces", "stress", "md_step"):
        fast_arr, slow_arr = getattr(fast, field), getattr(slow, field)
        assert fast_arr is not None, field
        assert slow_arr is not None, field
        np.testing.assert_array_equal(fast_arr, slow_arr, err_msg=field)


def test_from_extxyz_real_format_stress_energy(tmp_path: Path) -> None:
    """A hand-written DynaMat-style file (3x3 asymmetric stress, free_energy alongside
    energy, leading-whitespace atom lines, extra per-atom columns) matches ASE; the
    asymmetric stress pins the Voigt index order (ASE reshapes the 9-vector
    Fortran-order, then takes [00, 11, 22, 12, 02, 01]).
    """
    lattice = '"4.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 4.0"'
    props = "species:S:1:pos:R:3:forces:R:3:momenta:R:3:tags:I:1"
    content = (
        f"2\nLattice={lattice} Properties={props} energy=-1.5 "
        'stress="1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0" free_energy=-1.6 pbc="T T T"\n'
        "Cu       0.10000       0.20000       0.30000      0.01      0.02      0.03"
        "  1 1 1  7\n"
        "Cu       1.10000       1.20000       1.30000      0.11      0.12      0.13"
        "  2 2 2  8\n"
        f"2\nLattice={lattice} Properties={props} energy=-2.5 "
        'stress="9.0 8.0 7.0 6.0 5.0 4.0 3.0 2.0 1.0" free_energy=-2.6 pbc="T T T"\n'
        "Cu       0.40000       0.50000       0.60000      0.04      0.05      0.06"
        "  3 3 3  9\n"
        "Cu       1.40000       1.50000       1.60000      0.14      0.15      0.16"
        "  4 4 4  10\n"
    )
    path = tmp_path / "traj.extxyz"
    path.write_text(content)

    fast = Trajectory.from_extxyz(str(path))
    slow = Trajectory.from_ase(read_trajectory(str(path)))
    assert fast.stress is not None
    np.testing.assert_array_equal(fast.stress, slow.stress)
    np.testing.assert_array_equal(fast.stress[0], [1, 5, 9, 8, 7, 4])
    np.testing.assert_array_equal(fast.energy, [-1.5, -2.5])  # energy, not free_energy
    np.testing.assert_array_equal(fast.positions, slow.positions)
    np.testing.assert_array_equal(fast.forces, slow.forces)

    # a frame lacking a result the others carry fails loud instead of dropping it
    first_frame = content[: content.index("2\nLattice", 1)]
    path.write_text(first_frame + first_frame.replace(" energy=-1.5", ""))
    with pytest.raises(ValueError, match="frame 1 is missing a comment field"):
        Trajectory.from_extxyz(str(path))


@pytest.mark.parametrize(
    ("frame1_kwargs", "match"),
    [
        # a different number of atom lines fails the frame-stride check
        ({"count": 3, "n_lines": 3}, "not a multiple of frame stride"),
        ({"species": "Ag"}, "species differ"),  # same count, other composition
        ({"pbc": "T T F"}, "pbc differs"),
        ({"props": "species:S:1:pos:R:3:forces:R:3"}, "Properties header differs"),
        # count line disagrees with frame 0 while the frame stride still aligns
        ({"count": 5}, "count line says 5 atoms"),
        # md_step on frame 1 only: partial md_step fails loud like from_ase
        ({"extra": " md_step=10"}, "frame 0 is missing a comment field"),
    ],
    ids=["varying_count", "species", "pbc", "properties", "count_line", "md_step"],
)
def test_from_extxyz_rejects_inconsistent_metadata(
    tmp_path: Path, frame1_kwargs: dict[str, Any], match: str
) -> None:
    """Frames inconsistent with frame 0 fail loudly instead of mixing systems."""
    path = tmp_path / "traj.extxyz"
    path.write_text(extxyz_frame() + extxyz_frame(**frame1_kwargs))
    with pytest.raises(ValueError, match=match):
        Trajectory.from_extxyz(str(path))
