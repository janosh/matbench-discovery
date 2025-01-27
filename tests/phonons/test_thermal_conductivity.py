"""Tests for thermal conductivity calculation module."""

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk
from ase.calculators.calculator import Calculator
from ase.calculators.emt import EMT
from phono3py.api_phono3py import Phono3py
from phonopy.structure.atoms import PhonopyAtoms
from pymatviz.enums import Key

from matbench_discovery.enums import MbdKey
from matbench_discovery.phonons import thermal_conductivity as ltc
from matbench_discovery.phonons.thermal_conductivity import calculate_fc2_set

NP_RNG = np.random.default_rng(seed=0)


@pytest.fixture
def test_atoms() -> Atoms:
    """Create a simple cubic test structure using Al (supported by EMT)."""
    atoms = bulk("Al", "fcc", a=4.05)
    atoms.info["fc2_supercell"] = 2 * np.eye(3)
    atoms.info["fc3_supercell"] = np.eye(3)
    atoms.info["q_mesh"] = [2, 2, 2]
    return atoms


@pytest.fixture
def test_ph3(test_atoms: Atoms) -> Phono3py:
    """Create a test Phono3py object."""
    return ltc.init_phono3py(test_atoms)


@pytest.fixture
def test_calculator() -> EMT:
    """Create a simple EMT calculator for testing."""
    return EMT()


def test_init_phono3py(test_atoms: Atoms) -> None:
    """Test initialization of Phono3py object."""
    ph3 = ltc.init_phono3py(test_atoms)
    assert isinstance(ph3, Phono3py)
    assert list(ph3.mesh_numbers) == [2, 2, 2]
    assert ph3.supercell_matrix.tolist() == np.eye(3).tolist()


def test_init_phono3py_missing_info(test_atoms: Atoms) -> None:
    """Test initialization with missing metadata."""
    test_atoms.info.pop("fc2_supercell")
    with pytest.raises(ValueError, match="fc2_supercell.*not found"):
        ltc.init_phono3py(test_atoms)


def test_calculate_fc2_set(test_ph3: Phono3py, test_calculator: EMT) -> None:
    """Test calculation of 2nd order force constants."""
    forces = ltc.calculate_fc2_set(test_ph3, test_calculator, pbar_kwargs={})
    assert isinstance(forces, np.ndarray)
    assert forces.ndim == 3  # (n_displacements, n_atoms, 3)
    assert forces.shape[-1] == 3


def test_calculate_fc3_set(test_ph3: Phono3py, test_calculator: EMT) -> None:
    """Test calculation of 3rd order force constants."""
    forces = ltc.calculate_fc3_set(test_ph3, test_calculator, pbar_kwargs={})
    assert isinstance(forces, np.ndarray)
    assert forces.ndim == 3
    assert forces.shape[-1] == 3


def test_get_fc2_and_freqs(test_ph3: Phono3py, test_calculator: EMT) -> None:
    """Test getting force constants and frequencies."""
    ph3, fc2_set, freqs = ltc.get_fc2_and_freqs(
        test_ph3, test_calculator, pbar_kwargs={}
    )
    assert isinstance(ph3, Phono3py)
    assert isinstance(fc2_set, np.ndarray)
    assert isinstance(freqs, np.ndarray)
    assert freqs.ndim == 2  # (n_qpoints, n_bands)
    n_bz_grid, n_bands = 15, 3
    assert freqs.shape == (n_bz_grid, n_bands)


def test_load_force_sets(test_ph3: Phono3py) -> None:
    """Test loading pre-computed force sets.

    We need to match the exact force set structure that Phono3py expects:
    - One force set per displacement
    - Forces must be in the correct shape for the supercell
    - Forces must be consistent with the symmetry of the crystal
    """
    # Get the actual number of atoms in the supercell
    n_atoms_fc2 = len(test_ph3.phonon_supercell)
    n_atoms_fc3 = len(test_ph3.supercell)
    # Create arrays of zeros for fc2 and fc3 force sets
    n_fc2_disp = len(test_ph3.phonon_supercells_with_displacements)
    n_fc3_disp = len(test_ph3.supercells_with_displacements)
    # Initialize force arrays
    fc2_set = np.zeros((n_fc2_disp, n_atoms_fc2, 3))
    fc3_set = np.zeros((n_fc3_disp, n_atoms_fc3, 3))

    # Create masks for non-None displacements
    fc2_mask = [sc is not None for sc in test_ph3.phonon_supercells_with_displacements]
    fc3_mask = [sc is not None for sc in test_ph3.supercells_with_displacements]

    # Generate random forces for displaced atoms
    fc2_set[fc2_mask, 0] = NP_RNG.random((sum(fc2_mask), 3)) * 0.1
    fc3_set[fc3_mask, 0] = NP_RNG.random((sum(fc3_mask), 3)) * 0.1

    # Set the forces directly on the Phono3py object
    test_ph3.phonon_forces = fc2_set
    test_ph3.forces = fc3_set

    # Now load the force sets
    ph3 = ltc.load_force_sets(test_ph3, fc2_set, fc3_set)
    assert isinstance(ph3, Phono3py)
    np.testing.assert_array_equal(ph3.phonon_forces, fc2_set)
    np.testing.assert_array_equal(ph3.forces, fc3_set)


@pytest.mark.parametrize("temperature", [300, 500, 1000])
def test_calculate_conductivity(
    test_ph3: Phono3py, test_calculator: EMT, temperature: float
) -> None:
    """Test thermal conductivity calculation."""
    # First need to compute force constants
    test_ph3, fc2_set, _ = ltc.get_fc2_and_freqs(
        test_ph3, test_calculator, pbar_kwargs={}
    )

    fc3_set = ltc.calculate_fc3_set(test_ph3, calculator=test_calculator)

    test_ph3.produce_fc3(symmetrize_fc3r=True)

    test_ph3 = ltc.load_force_sets(test_ph3, fc2_set, fc3_set)

    ph3, kappa_dict, kappa = ltc.calculate_conductivity(test_ph3, [temperature])

    assert isinstance(ph3, Phono3py)
    assert isinstance(kappa_dict, dict)
    required_keys = {
        MbdKey.kappa_tot_rta,
        MbdKey.kappa_p_rta,
        MbdKey.kappa_c,
        Key.mode_weights,
        Key.q_points,
        Key.ph_freqs,
        MbdKey.mode_kappa_tot,
    }
    assert set(kappa_dict) >= required_keys
    assert all(isinstance(val, np.ndarray) for val in kappa_dict.values())


def test_calculate_mode_kappa_tot() -> None:
    """Test calculation of total mode kappa."""
    # Create dummy arrays with correct shapes
    mode_kappa_p_rta = NP_RNG.random((2, 3, 3, 3))  # (T, q-points, bands, xyz)
    mode_kappa_coherence = NP_RNG.random(
        (2, 3, 3, 3, 3)
    )  # (T, q-points, bands, bands, xyz)
    heat_capacity = NP_RNG.random((2, 3, 3))  # (T, q-points, bands)

    result = ltc.calculate_mode_kappa_tot(
        mode_kappa_p_rta, mode_kappa_coherence, heat_capacity
    )
    assert isinstance(result, np.ndarray)
    assert result.shape == mode_kappa_p_rta.shape


class MockCalculator(Calculator):
    """Mock calculator that returns predefined forces."""

    def __init__(self, forces: np.ndarray) -> None:
        super().__init__()
        self.forces = forces

    def get_forces(self, _atoms: Atoms) -> np.ndarray:
        """Return predefined forces."""
        return self.forces


def test_calculate_fc2_set_forces() -> None:
    """Test that calculate_fc2_set correctly handles force calculations."""
    # Create a simple cubic structure using PhonopyAtoms
    atoms = PhonopyAtoms(
        symbols=["Si"] * 2,
        scaled_positions=[[0, 0, 0], [0.25, 0.25, 0.25]],
        cell=[[1, 0, 0], [0, 1, 0], [0, 1, 1]],
    )

    # Initialize Phono3py object with the test structure
    ph3 = Phono3py(atoms, supercell_matrix=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    # Generate the displacements and create dataset
    ph3.generate_displacements(distance=0.03)
    ph3.generate_fc2_displacements()

    # Create mock forces that would be physically reasonable
    # Forces sum to zero (Newton's 3rd law)
    mock_forces = np.hstack((0.1 * np.eye(3), -0.1 * np.eye(3)))

    # Test with mock calculator that returns our predefined forces
    calc = MockCalculator(mock_forces[0])  # Start with first set of forces
    force_set = calculate_fc2_set(ph3, calc, pbar_kwargs={"disable": True})

    # Test shape of returned force set
    expected_shape = (len(ph3.phonon_supercells_with_displacements), 6)
    assert force_set.shape == expected_shape, f"{force_set.shape=} != {expected_shape=}"

    # Test that forces sum to zero for each displacement (conservation of momentum)
    for forces in force_set:
        assert np.allclose(forces.sum(axis=0), [0, 0, 0], atol=1e-10), (
            "Forces don't sum to zero"
        )

    # Test that using np.ones instead of actual forces would break conservation laws
    wrong_forces = np.ones_like(mock_forces[0])
    calc_wrong = MockCalculator(wrong_forces)
    wrong_force_set = calculate_fc2_set(ph3, calc_wrong, pbar_kwargs={"disable": True})

    # This should not sum to zero as all forces are 1
    assert not np.allclose(wrong_force_set.sum(axis=1), 0), (
        "Incorrect forces (all ones) incorrectly sum to zero"
    )


def test_calculate_fc2_set_null_supercell() -> None:
    """Test that calculate_fc2_set handles null supercells correctly."""
    atoms = PhonopyAtoms(
        symbols=["Si"] * 2,
        scaled_positions=[[0, 0, 0], [0.25, 0.25, 0.25]],
        cell=[[1, 0, 0], [0, 1, 0], [0, 1, 1]],
    )
    ph3 = Phono3py(atoms, supercell_matrix=np.eye(3))

    # Generate displacements and create dataset
    ph3.generate_displacements(distance=0.03)
    ph3.generate_fc2_displacements()

    calc = MockCalculator(np.zeros((2, 3)))
    force_set = calculate_fc2_set(ph3, calc, pbar_kwargs={"disable": True})

    # Test that ones are used for null supercell case
    np.testing.assert_allclose(force_set, np.zeros((2, len(ph3.phonon_supercell), 3)))

    # Verify shape is correct even with null supercell
    expected_shape = (2, len(ph3.phonon_supercell), 3)
    assert force_set.shape == expected_shape, f"{force_set.shape=} != {expected_shape=}"
