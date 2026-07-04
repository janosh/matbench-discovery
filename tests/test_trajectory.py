"""Tests for the dense-array Trajectory container and its HDF5 / ASE interop."""

from pathlib import Path

import numpy as np
import pytest
from ase import Atoms

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


def make_atoms(symbols: str = "Cu2", *, pbc: bool = True) -> Atoms:
    """Minimal ASE Atoms (positions at the origin) in a 4 A cubic cell."""
    return Atoms(symbols, cell=np.eye(3) * 4, pbc=pbc)


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
    traj.write_hdf5(path)

    full = Trajectory.read_hdf5(path)
    np.testing.assert_array_equal(full.positions, traj.positions)
    np.testing.assert_array_equal(full.cell, traj.cell)
    assert (full.atomic_numbers == traj.atomic_numbers).all()
    assert (full.pbc == traj.pbc).all()
    for field in ("energy", "forces", "stress", "md_step"):
        if results:
            np.testing.assert_array_equal(getattr(full, field), getattr(traj, field))
        else:
            assert getattr(full, field) is None

    strided = Trajectory.read_hdf5(path, frames=slice(0, 20, 4))
    assert strided.n_frames == 5
    np.testing.assert_array_equal(strided.positions, traj.positions[0:20:4])


def test_read_hdf5_schema_mismatch(tmp_path: Path) -> None:
    """A trajectory file with an unexpected schema fails closed."""
    import h5py

    path = str(tmp_path / "traj.h5")
    make_traj().write_hdf5(path)
    with h5py.File(path, "r+") as file:
        file.attrs["schema"] = TRAJECTORY_SCHEMA + 1
    with pytest.raises(ValueError, match="trajectory schema"):
        Trajectory.read_hdf5(path)


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
    ],
    ids=["empty", "ragged_count", "species", "pbc"],
)
def test_from_ase_rejects_bad_input(frames: list[Atoms], match: str) -> None:
    """from_ase rejects empty input and frames inconsistent with frame 0."""
    with pytest.raises(ValueError, match=match):
        Trajectory.from_ase(frames)


def test_from_ase_reraises_genuine_backend_error() -> None:
    """from_ase returns None for unavailable properties but re-raises real failures."""
    from ase.calculators.calculator import Calculator, all_changes

    class BrokenStress(Calculator):
        implemented_properties = ("stress",)

        def calculate(
            self,
            atoms: Atoms | None = None,
            properties: list[str] | None = None,
            system_changes: list[str] = all_changes,
        ) -> None:
            super().calculate(atoms, properties, system_changes)
            raise RuntimeError("simulated stress backend failure")

    atoms = make_atoms()
    atoms.calc = BrokenStress()
    with pytest.raises(RuntimeError, match="stress backend failure"):
        Trajectory.from_ase([atoms])


def test_from_ase_rejects_mixed_property_availability() -> None:
    """A property present on some frames but not others fails loud (not silent None)."""
    from ase.calculators.singlepoint import SinglePointCalculator

    with_energy = make_atoms()
    with_energy.calc = SinglePointCalculator(with_energy, energy=-1.0)
    with pytest.raises(ValueError, match="missing a property"):
        Trajectory.from_ase([with_energy, make_atoms()])

    # md_step on some frames but not others must also fail loud (not silently drop it)
    with_step = make_atoms()
    with_step.info["md_step"] = 0
    with pytest.raises(ValueError, match="missing md_step"):
        Trajectory.from_ase([with_step, make_atoms()])


@pytest.mark.parametrize("suffix", ["", ".xz", ".gz"])
def test_from_extxyz_matches_ase_parser(tmp_path: Path, suffix: str) -> None:
    """The fast bulk reader reproduces the ASE extxyz reader, which also round-trips
    the original trajectory within extxyz text precision, including compression.
    """
    import ase.io

    from matbench_discovery.md import read_trajectory

    traj = make_traj(n_frames=5, n_atoms=4)
    assert traj.energy is not None
    # non-orthogonal cell so the parity check catches Lattice row/column-order or
    # transpose bugs that a diagonal cell would hide
    traj.cell[:] = [[5.0, 0.5, 0.3], [0.0, 5.2, 0.4], [0.0, 0.0, 4.7]]
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
    energy, leading-whitespace atom lines) matches ASE; the asymmetric stress pins the
    Voigt index order (ASE reshapes the 9-vector Fortran-order, then takes
    [00, 11, 22, 12, 02, 01]).
    """
    from matbench_discovery.md import read_trajectory

    lattice = '"4.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 4.0"'
    props = "species:S:1:pos:R:3:forces:R:3"
    content = (
        f"2\nLattice={lattice} Properties={props} energy=-1.5 "
        'stress="1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0" free_energy=-1.6 pbc="T T T"\n'
        "Cu       0.10000       0.20000       0.30000      0.01      0.02      0.03\n"
        "Cu       1.10000       1.20000       1.30000      0.11      0.12      0.13\n"
        f"2\nLattice={lattice} Properties={props} energy=-2.5 "
        'stress="9.0 8.0 7.0 6.0 5.0 4.0 3.0 2.0 1.0" free_energy=-2.6 pbc="T T T"\n'
        "Cu       0.40000       0.50000       0.60000      0.04      0.05      0.06\n"
        "Cu       1.40000       1.50000       1.60000      0.14      0.15      0.16\n"
    )
    path = str(tmp_path / "traj.extxyz")
    with open(path, "w") as file:
        file.write(content)

    fast = Trajectory.from_extxyz(path)
    slow = Trajectory.from_ase(read_trajectory(path))
    assert fast.stress is not None
    np.testing.assert_array_equal(fast.stress, slow.stress)
    np.testing.assert_array_equal(fast.stress[0], [1, 5, 9, 8, 7, 4])
    np.testing.assert_array_equal(fast.energy, [-1.5, -2.5])  # energy, not free_energy
    np.testing.assert_array_equal(fast.positions, slow.positions)
    np.testing.assert_array_equal(fast.forces, slow.forces)


@pytest.mark.parametrize(
    ("frame1", "match"),
    [
        # frame 1 has a different number of atom lines -> fails the frame-stride check
        (
            '3\nLattice="4 0 0 0 4 0 0 0 4" '
            'Properties=species:S:1:pos:R:3 pbc="T T T"\n'
            "Cu 0 0 0\nCu 1 1 1\nCu 2 2 2\n",
            "not a multiple of frame stride",
        ),
        # same atom count but different species composition
        (
            '2\nLattice="4 0 0 0 4 0 0 0 4" '
            'Properties=species:S:1:pos:R:3 pbc="T T T"\n'
            "Ag 0 0 0\nAg 1 1 1\n",
            "species differ",
        ),
        # different periodic boundary conditions
        (
            '2\nLattice="4 0 0 0 4 0 0 0 4" '
            'Properties=species:S:1:pos:R:3 pbc="T T F"\n'
            "Cu 0 0 0\nCu 1 1 1\n",
            "pbc differs",
        ),
        # different Properties header (column layout)
        (
            '2\nLattice="4 0 0 0 4 0 0 0 4" '
            'Properties=species:S:1:pos:R:3:forces:R:3 pbc="T T T"\n'
            "Cu 0 0 0\nCu 1 1 1\n",
            "Properties header differs",
        ),
        # count line value disagrees with frame 0 while the frame stride still aligns
        (
            '5\nLattice="4 0 0 0 4 0 0 0 4" '
            'Properties=species:S:1:pos:R:3 pbc="T T T"\n'
            "Cu 0 0 0\nCu 1 1 1\n",
            "count line says 5 atoms",
        ),
    ],
    ids=["varying_count", "species", "pbc", "properties", "count_line"],
)
def test_from_extxyz_rejects_inconsistent_metadata(
    tmp_path: Path, frame1: str, match: str
) -> None:
    """Frames differing from frame 0 in atom count, species, Properties header or pbc
    must fail loudly rather than silently mixing systems into the stacked arrays.
    """
    frame0 = (
        '2\nLattice="4 0 0 0 4 0 0 0 4" '
        'Properties=species:S:1:pos:R:3 pbc="T T T"\n'
        "Cu 0 0 0\nCu 1 1 1\n"
    )
    path = str(tmp_path / "traj.extxyz")
    with open(path, "w") as file:
        file.write(frame0 + frame1)
    with pytest.raises(ValueError, match=match):
        Trajectory.from_extxyz(path)
