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
    """len/n_frames/n_atoms and frame slicing return consistent sub-trajectories."""
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


def test_from_ase_missing_results_and_bad_input() -> None:
    """from_ase yields None for absent results and rejects empty/ragged input."""
    bare = [Atoms("Cu2", positions=np.zeros((2, 3)), cell=np.eye(3) * 4, pbc=True)]
    traj = Trajectory.from_ase(bare)
    assert traj.energy is None
    assert traj.forces is None
    assert traj.stress is None
    assert traj.md_step is None

    with pytest.raises(ValueError, match="zero frames"):
        Trajectory.from_ase([])
    ragged = [
        Atoms("Cu2", positions=np.zeros((2, 3)), cell=np.eye(3) * 4, pbc=True),
        Atoms("Cu3", positions=np.zeros((3, 3)), cell=np.eye(3) * 4, pbc=True),
    ]
    with pytest.raises(ValueError, match="Inconsistent atom counts"):
        Trajectory.from_ase(ragged)


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

    atoms = Atoms("Cu2", positions=np.zeros((2, 3)), cell=np.eye(3) * 4, pbc=True)
    atoms.calc = BrokenStress()
    with pytest.raises(RuntimeError, match="stress backend failure"):
        Trajectory.from_ase([atoms])


def test_from_ase_parses_extxyz(tmp_path: Path) -> None:
    """from_ase rebuilds a Trajectory from ASE-parsed extxyz (the converter path)."""
    import ase.io

    from matbench_discovery.md import read_trajectory

    traj = make_traj(n_frames=4, n_atoms=3)
    assert traj.energy is not None
    path = str(tmp_path / "traj.extxyz")
    ase.io.write(path, traj.to_ase())

    parsed = Trajectory.from_ase(read_trajectory(path))
    assert parsed.energy is not None
    assert parsed.n_frames == 4
    # extxyz is a text format, so positions/energy round-trip to its print precision
    np.testing.assert_allclose(parsed.positions, traj.positions, rtol=0, atol=1e-6)
    np.testing.assert_allclose(parsed.energy, traj.energy, rtol=0, atol=1e-6)


def test_from_ase_rejects_inconsistent_species_or_pbc() -> None:
    """from_ase requires each frame to match frame 0 species and pbc, not just count."""
    cu = Atoms("Cu2", positions=np.zeros((2, 3)), cell=np.eye(3) * 4, pbc=True)
    he = Atoms("He2", positions=np.zeros((2, 3)), cell=np.eye(3) * 4, pbc=True)
    with pytest.raises(ValueError, match="species differ"):
        Trajectory.from_ase([cu, he])
    mixed_pbc = Atoms("Cu2", positions=np.zeros((2, 3)), cell=np.eye(3) * 4, pbc=False)
    with pytest.raises(ValueError, match="pbc differs"):
        Trajectory.from_ase([cu, mixed_pbc])


def test_from_ase_rejects_mixed_property_availability() -> None:
    """A property present on some frames but not others fails loud (not silent None)."""
    from ase.calculators.singlepoint import SinglePointCalculator

    with_energy = Atoms("Cu2", positions=np.zeros((2, 3)), cell=np.eye(3) * 4, pbc=True)
    with_energy.calc = SinglePointCalculator(with_energy, energy=-1.0)
    without = Atoms("Cu2", positions=np.zeros((2, 3)), cell=np.eye(3) * 4, pbc=True)
    with pytest.raises(ValueError, match="missing a property"):
        Trajectory.from_ase([with_energy, without])


def test_negative_and_out_of_range_indexing() -> None:
    """Integer indexing supports negative indices and bounds-checks."""
    traj = make_traj(n_frames=5, n_atoms=2)
    np.testing.assert_array_equal(traj[-1].positions[0], traj.positions[4])
    np.testing.assert_array_equal(traj.frame_as_atoms(-1).positions, traj.positions[4])
    with pytest.raises(IndexError):
        _ = traj[5]
    with pytest.raises(IndexError):
        _ = traj[-6]


def test_to_ase_dropped_calc_on_missing_results() -> None:
    """to_ase omits a calculator when the trajectory carries no results."""
    atoms = make_traj(results=False).to_ase()[0]
    assert atoms.calc is None
