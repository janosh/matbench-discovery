"""Tests for thermal conductivity calculation module."""

from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk
from ase.calculators.calculator import Calculator
from ase.calculators.emt import EMT
from phono3py.api_phono3py import Phono3py
from phono3py.conductivity.calculators import RTACalculator
from phonopy.structure.atoms import PhonopyAtoms
from pymatviz.enums import Key

from matbench_discovery.enums import MbdKey
from matbench_discovery.phonons import thermal_conductivity as ltc
from matbench_discovery.phonons.thermal_conductivity import (
    calc_kappa_for_structure,
    calculate_fc2_set,
)

NP_RNG = np.random.default_rng(seed=0)


@pytest.fixture
def test_atoms() -> Atoms:
    """Create a simple FCC test structure using Cu (supported by EMT).

    Cu is used instead of Al because Al is mono-isotopic and triggers a
    phonopy bug where single-isotope entries lack a trailing comma, causing
    ``((27, 26.98, 1.0))`` to be parsed as a flat tuple instead of a nested one.
    """
    atoms = bulk("Cu", "fcc", a=3.615)
    atoms.info["fc2_supercell"] = 2 * np.eye(3)
    atoms.info["fc3_supercell"] = np.eye(3)
    atoms.info["q_point_mesh"] = (2, 2, 2)
    return atoms


@pytest.fixture
def test_ph3(test_atoms: Atoms) -> Phono3py:
    """Create a test Phono3py object."""
    return ltc.init_phono3py(
        atoms=test_atoms,
        fc2_supercell=test_atoms.info["fc2_supercell"],
        fc3_supercell=test_atoms.info["fc3_supercell"],
        q_point_mesh=test_atoms.info["q_point_mesh"],
    )


@pytest.fixture
def test_calculator() -> EMT:
    """Create a simple EMT calculator for testing."""
    return EMT()


def test_init_phono3py(test_atoms: Atoms) -> None:
    """Test initialization of Phono3py object."""
    fc2_supercell = test_atoms.info["fc2_supercell"]
    fc3_supercell = test_atoms.info["fc3_supercell"]
    q_point_mesh = test_atoms.info["q_point_mesh"]
    ph3 = ltc.init_phono3py(
        test_atoms,
        fc2_supercell=fc2_supercell,
        fc3_supercell=fc3_supercell,
        q_point_mesh=q_point_mesh,
    )
    assert isinstance(ph3, Phono3py)
    assert ph3.mesh_numbers is not None
    assert ph3.mesh_numbers.tolist() == [2, 2, 2]
    assert ph3.supercell_matrix.tolist() == np.eye(3).tolist()
    # Check that both supercells were created correctly
    assert ph3.phonon_supercell_matrix is not None
    assert np.allclose(ph3.phonon_supercell_matrix, fc2_supercell)
    assert np.allclose(ph3.supercell_matrix, fc3_supercell)
    # Verify supercell sizes
    assert len(ph3.phonon_supercell) == len(test_atoms) * 8  # 2x2x2 supercell
    assert len(ph3.supercell) == len(test_atoms)  # 1x1x1 supercell

    disp_magnitudes = np.linalg.norm(ph3.displacements, axis=-1)
    # Assert that default displacement distance is 0.01
    np.testing.assert_allclose(disp_magnitudes, 0.01)


def test_init_phono3py_custom_mesh_and_displacement(test_atoms: Atoms) -> None:
    """Test initialization with custom mesh and displacement distance."""
    fc2_supercell = test_atoms.info["fc2_supercell"]
    fc3_supercell = test_atoms.info["fc3_supercell"]
    q_point_mesh = test_atoms.info["q_point_mesh"]

    # Test custom mesh
    custom_mesh = (4, 4, 4)
    ph3 = ltc.init_phono3py(
        test_atoms,
        fc2_supercell=fc2_supercell,
        fc3_supercell=fc3_supercell,
        q_point_mesh=custom_mesh,
    )
    assert ph3.mesh_numbers is not None
    assert tuple(ph3.mesh_numbers.tolist()) == custom_mesh
    # Phono3py automatically initializes mesh when setting mesh_numbers
    assert (bz_grid := ph3.grid) is not None
    assert len(bz_grid.addresses) == 89

    custom_displacement = 0.05
    ph3 = ltc.init_phono3py(
        test_atoms,
        fc2_supercell=fc2_supercell,
        fc3_supercell=fc3_supercell,
        q_point_mesh=q_point_mesh,
        displacement_distance=custom_displacement,
    )
    # Check that displacements were generated with the custom distance
    disp_magnitudes = np.linalg.norm(ph3.displacements, axis=-1)
    non_zero_disps = disp_magnitudes[disp_magnitudes > 0]
    np.testing.assert_allclose(non_zero_disps, custom_displacement)


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
        MbdKey.mode_kappa_tot_rta,
    }
    assert set(kappa_dict) >= required_keys
    assert all(isinstance(val, np.ndarray) for val in kappa_dict.values())
    assert isinstance(kappa, RTACalculator)
    assert (freqs := kappa.frequencies) is not None
    assert freqs.shape == (3, 3)
    assert kappa.grid_points.shape == (3,)
    assert (gamma := kappa.gamma) is not None
    assert gamma.shape == (1, 1, 3, 3)


def test_calculate_mode_kappa_tot() -> None:
    """Test calculation of total mode kappa."""
    # Create dummy arrays with correct shapes
    mode_kappa_p_rta = NP_RNG.random((2, 3, 3, 3))  # (T, q-points, bands, xyz)
    mode_kappa_coherence = NP_RNG.random(
        (2, 3, 3, 3, 3)
    )  # (T, q-points, bands, bands, xyz)
    heat_capacity = NP_RNG.random((2, 3, 3))  # (T, q-points, bands)

    result = ltc.calc_mode_kappa_tot(
        mode_kappa_p_rta, mode_kappa_coherence, heat_capacity
    )
    assert isinstance(result, np.ndarray)
    assert result.shape == mode_kappa_p_rta.shape


class MockCalculator(Calculator):
    """Mock calculator that returns predefined forces."""

    def __init__(self, forces: np.ndarray) -> None:
        """Store the forces returned by get_forces()."""
        super().__init__()
        self.forces = forces

    def get_forces(self, atoms: Atoms | None = None) -> np.ndarray:  # noqa: ARG002
        """Return predefined forces."""
        return self.forces


def make_si2_phonopy_atoms() -> PhonopyAtoms:
    """Create a minimal Si2 diamond-like cell as PhonopyAtoms."""
    return PhonopyAtoms(
        symbols=["Si", "Si"],
        scaled_positions=[[0, 0, 0], [0.25, 0.25, 0.25]],
        cell=[[1, 0, 0], [0, 1, 0], [0, 1, 1]],
    )


def build_ph3_with_fc2(
    atoms: PhonopyAtoms,
    fc2_matrix: np.ndarray,
    *,
    fc3_matrix: np.ndarray | None = None,
    distance: float = 0.02,
) -> Phono3py:
    """Construct Phono3py with a given phonon supercell and make fc2 displacements."""
    ph3_local = Phono3py(
        atoms,
        supercell_matrix=fc3_matrix if fc3_matrix is not None else np.eye(3, dtype=int),
        phonon_supercell_matrix=fc2_matrix,
    )
    ph3_local.generate_displacements(distance=distance)
    ph3_local.generate_fc2_displacements()
    return ph3_local


def test_calculate_fc2_set_forces() -> None:
    """Test that calculate_fc2_set correctly handles force calculations."""
    atoms = make_si2_phonopy_atoms()

    # Initialize Phono3py object with the test structure
    ph3 = build_ph3_with_fc2(atoms, np.eye(3, dtype=int), distance=0.03)

    # Create mock forces that would be physically reasonable
    # Forces sum to zero (Newton's 3rd law)
    mock_forces = np.hstack((0.1 * np.eye(3), -0.1 * np.eye(3)))

    # Test with mock calculator that returns our predefined forces
    # Provide forces in shape (n_atoms, 3)
    calc = MockCalculator(mock_forces[0].reshape(2, 3))
    force_set = calculate_fc2_set(ph3, calc, pbar_kwargs={"disable": True})

    # Test shape of returned force set
    n_displacements = len(ph3.phonon_supercells_with_displacements)
    n_atoms = len(ph3.phonon_supercell)
    expected_shape = (n_displacements, n_atoms, 3)
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
    atoms = make_si2_phonopy_atoms()
    ph3 = build_ph3_with_fc2(atoms, np.eye(3, dtype=int), distance=0.03)

    calc = MockCalculator(np.zeros((2, 3)))
    force_set = calculate_fc2_set(ph3, calc, pbar_kwargs={"disable": True})

    # Test that zeros are used for null supercell case
    all_zeros = np.zeros(
        (len(ph3.phonon_supercells_with_displacements), len(ph3.phonon_supercell), 3)
    )
    np.testing.assert_allclose(force_set, all_zeros)

    # Verify shape is correct even with null supercell
    expected_shape = (
        len(ph3.phonon_supercells_with_displacements),
        len(ph3.phonon_supercell),
        3,
    )
    assert force_set.shape == expected_shape, f"{force_set.shape=} != {expected_shape=}"


@pytest.mark.parametrize(
    "fc2_matrix, expected_det",
    [
        (np.eye(3, dtype=int), 1),
        (2 * np.eye(3, dtype=int), 8),
        (np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]], dtype=int), 1),
        (np.array([[2, 1, 0], [0, 1, 0], [0, 0, 1]], dtype=int), 2),
    ],
)
def test_calculate_fc2_set_with_various_supercells(
    fc2_matrix: np.ndarray, expected_det: int
) -> None:
    """Check fc2 forces shape across varied phonon supercells."""
    atoms = make_si2_phonopy_atoms()
    ph3 = build_ph3_with_fc2(atoms, fc2_matrix)

    # Calculator returns zeros with correct per-supercell shape
    calc = MockCalculator(np.zeros((len(ph3.phonon_supercell), 3)))
    force_set = calculate_fc2_set(ph3, calc, pbar_kwargs={"disable": True})

    # Expected atoms in phonon supercell is n_atoms * det(fc2_matrix)
    n_atoms = len(atoms)
    assert len(ph3.phonon_supercell) == n_atoms * expected_det

    expected_shape = (
        len(ph3.phonon_supercells_with_displacements),
        len(ph3.phonon_supercell),
        3,
    )
    assert force_set.shape == expected_shape


def test_calculate_fc2_set_requires_phonon_supercell() -> None:
    """Calling fc2 forces without phonon supercell must raise."""
    atoms = make_si2_phonopy_atoms()
    ph3 = Phono3py(atoms, supercell_matrix=np.eye(3, dtype=int))
    ph3.generate_displacements(distance=0.02)
    ph3.generate_fc2_displacements()
    calc = MockCalculator(np.zeros((len(ph3.supercell), 3)))
    with pytest.raises(RuntimeError):
        _ = calculate_fc2_set(ph3, calc, pbar_kwargs={"disable": True})


# === calc_kappa_for_structure: per-structure relax -> FC -> kappa orchestrator ===


@pytest.fixture
def mock_atoms() -> Atoms:
    """Create a simple H2 atoms object for testing."""
    atoms = Atoms("H2", positions=[[0, 0, 0], [0.74, 0, 0]], cell=[10, 10, 10])
    atoms.info[Key.mat_id] = "test-mat-id"
    atoms.info["fc2_supercell"] = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    atoms.info["fc3_supercell"] = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    atoms.info["q_point_mesh"] = [2, 2, 2]
    return atoms


@pytest.fixture
def mock_calculator() -> MagicMock:
    """Create a mock calculator with standard returns."""
    calc = MagicMock()
    calc.get_potential_energy.return_value = -1.0
    calc.get_forces.return_value = np.zeros((2, 3))
    calc.get_stress.return_value = np.zeros(6)
    return calc


@pytest.fixture
def common_kappa_kwargs(tmp_path_factory: pytest.TempPathFactory) -> dict:
    """Common kwargs for calc_kappa_for_structure calls."""
    return {
        "displacement_distance": 0.01,
        "temperatures": [300],
        "ase_optimizer": "FIRE",
        "max_steps": 10,
        "force_max": 0.01,
        "symprec": 1e-5,
        "enforce_relax_symm": False,
        "save_forces": False,
        "out_dir": str(tmp_path_factory.mktemp("test_kappa")),
        "task_id": 0,
    }


@pytest.mark.parametrize(
    ("init_spg", "relaxed_spg", "expected_broken"),
    [
        (139, 1, True),
        (1, 139, True),
        (139, 139, False),
        (225, 225, False),
        (16, 2, True),
    ],
)
def test_broken_symmetry_detection(
    mock_atoms: Atoms,
    mock_calculator: MagicMock,
    common_kappa_kwargs: dict,
    init_spg: int,
    relaxed_spg: int,
    expected_broken: bool,
) -> None:
    """Test symmetry breaking detection for various space group transitions."""
    freqs_2d = np.array([[0.001, 0.002, 0.003, 1.0, 2.0, 3.0]])
    with (
        patch("matbench_discovery.phonons.thermal_conductivity.MoyoDataset") as m_moyo,
        patch("matbench_discovery.phonons.thermal_conductivity.MoyoAdapter"),
        patch("matbench_discovery.phonons.thermal_conductivity.init_phono3py"),
        patch(
            "matbench_discovery.phonons.thermal_conductivity.get_fc2_and_freqs"
        ) as m_fc2,
        patch("matbench_discovery.phonons.thermal_conductivity.calculate_fc3_set"),
        patch("matbench_discovery.phonons.thermal_conductivity.calculate_conductivity"),
    ):
        m_moyo.side_effect = [
            MagicMock(number=init_spg),
            MagicMock(number=relaxed_spg),
        ]
        m_fc2.return_value = (MagicMock(), [], freqs_2d)

        mat_id, results, _ = calc_kappa_for_structure(
            atoms=mock_atoms, calculator=mock_calculator, **common_kappa_kwargs
        )

        assert mat_id == "test-mat-id"
        assert results["broken_symmetry"] is expected_broken
        assert results["relaxed_space_group_number"] == relaxed_spg


@pytest.mark.parametrize(
    ("has_imag", "broken", "allow_broken", "ignore_imag", "should_calc"),
    [
        # Standard cases: ignore_imaginary_freqs=False (MACE mode)
        (False, False, False, False, True),
        (False, True, False, False, False),
        (False, True, True, False, True),
        (True, False, False, False, False),
        (True, True, False, False, False),
        (True, True, True, False, False),
        # NequIP/Allegro mode: ignore_imaginary_freqs=True
        (False, False, False, True, True),
        (False, True, False, True, True),
        (True, False, False, True, True),
        (True, True, False, True, True),
    ],
)
def test_ltc_calculation_conditions(
    mock_atoms: Atoms,
    mock_calculator: MagicMock,
    common_kappa_kwargs: dict,
    has_imag: bool,
    broken: bool,
    allow_broken: bool,
    ignore_imag: bool,
    should_calc: bool,
) -> None:
    """Test LTC calculation proceeds only under correct conditions."""
    freqs = (
        np.array([[-1.0, 0.1, 0.2, 1.0, 2.0, 3.0]])
        if has_imag
        else np.array([[0.001, 0.002, 0.003, 1.0, 2.0, 3.0]])
    )

    with (
        patch("matbench_discovery.phonons.thermal_conductivity.MoyoDataset") as m_moyo,
        patch("matbench_discovery.phonons.thermal_conductivity.MoyoAdapter"),
        patch("matbench_discovery.phonons.thermal_conductivity.init_phono3py"),
        patch(
            "matbench_discovery.phonons.thermal_conductivity.get_fc2_and_freqs"
        ) as m_fc2,
        patch(
            "matbench_discovery.phonons.thermal_conductivity.calculate_fc3_set"
        ) as m_fc3,
    ):
        m_moyo.side_effect = [
            MagicMock(number=139),
            MagicMock(number=1 if broken else 139),
        ]
        m_fc2.return_value = (MagicMock(), [], freqs)
        m_fc3.return_value = []

        calc_kappa_for_structure(
            atoms=mock_atoms,
            calculator=mock_calculator,
            conductivity_broken_symm=allow_broken,
            ignore_imaginary_freqs=ignore_imag,
            **common_kappa_kwargs,
        )

        (m_fc3.assert_called_once if should_calc else m_fc3.assert_not_called)()


@pytest.mark.parametrize(
    ("max_steps", "should_check"),
    [(10, True), (0, False), (1, True), (1000, True)],
)
def test_symmetry_check_only_when_relaxing(
    mock_atoms: Atoms,
    mock_calculator: MagicMock,
    common_kappa_kwargs: dict,
    max_steps: int,
    should_check: bool,
) -> None:
    """Test that symmetry is only checked when relaxation occurs."""
    freqs_2d = np.array([[0.001, 0.002, 0.003, 1.0, 2.0, 3.0]])
    with (
        patch("matbench_discovery.phonons.thermal_conductivity.MoyoDataset") as m_moyo,
        patch("matbench_discovery.phonons.thermal_conductivity.MoyoAdapter"),
        patch("matbench_discovery.phonons.thermal_conductivity.init_phono3py"),
        patch(
            "matbench_discovery.phonons.thermal_conductivity.get_fc2_and_freqs"
        ) as m_fc2,
        patch("matbench_discovery.phonons.thermal_conductivity.calculate_fc3_set"),
        patch("matbench_discovery.phonons.thermal_conductivity.calculate_conductivity"),
    ):
        m_moyo.return_value = MagicMock(number=139)
        m_fc2.return_value = (MagicMock(), [], freqs_2d)

        kwargs = {k: v for k, v in common_kappa_kwargs.items() if k != "max_steps"}
        _, results, _ = calc_kappa_for_structure(
            atoms=mock_atoms,
            calculator=mock_calculator,
            max_steps=max_steps,
            **kwargs,
        )

        assert "broken_symmetry" in results
        if should_check:
            assert m_moyo.call_count == 2
            assert "relaxed_space_group_number" in results
        else:
            assert m_moyo.call_count == 1
            assert results["broken_symmetry"] is False


@pytest.mark.parametrize(
    ("symprec", "init_spg", "relaxed_spg", "expected"),
    [
        (1e-5, 139, 139, False),
        (1e-5, 139, 138, True),
        (1e-2, 139, 139, False),
        (1e-1, 225, 225, False),
    ],
)
def test_symprec_parameter(
    mock_atoms: Atoms,
    mock_calculator: MagicMock,
    common_kappa_kwargs: dict,
    symprec: float,
    init_spg: int,
    relaxed_spg: int,
    expected: bool,
) -> None:
    """Test that symprec parameter is correctly passed to MoyoDataset."""
    freqs_2d = np.array([[0.001, 0.002, 0.003, 1.0, 2.0, 3.0]])
    with (
        patch("matbench_discovery.phonons.thermal_conductivity.MoyoDataset") as m_moyo,
        patch("matbench_discovery.phonons.thermal_conductivity.MoyoAdapter"),
        patch("matbench_discovery.phonons.thermal_conductivity.init_phono3py"),
        patch(
            "matbench_discovery.phonons.thermal_conductivity.get_fc2_and_freqs"
        ) as m_fc2,
        patch("matbench_discovery.phonons.thermal_conductivity.calculate_fc3_set"),
        patch("matbench_discovery.phonons.thermal_conductivity.calculate_conductivity"),
    ):
        m_moyo.side_effect = [
            MagicMock(number=init_spg),
            MagicMock(number=relaxed_spg),
        ]
        m_fc2.return_value = (MagicMock(), [], freqs_2d)

        kwargs = {k: v for k, v in common_kappa_kwargs.items() if k != "symprec"}
        _, results, _ = calc_kappa_for_structure(
            atoms=mock_atoms,
            calculator=mock_calculator,
            symprec=symprec,
            **kwargs,
        )

        for call in m_moyo.call_args_list:
            assert call[1]["symprec"] == symprec
        assert results["broken_symmetry"] is expected
