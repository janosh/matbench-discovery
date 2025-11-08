"""Tests for calc_kappa.py."""

from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from ase import Atoms
from pymatviz.enums import Key


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
def common_kappa_kwargs() -> dict:
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
        "out_dir": "/tmp/test_kappa",
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
    from matbench_discovery.phonons.calc_kappa import calc_kappa_for_structure

    freqs_2d = np.array([[0.001, 0.002, 0.003, 1.0, 2.0, 3.0]])
    with (
        patch("matbench_discovery.phonons.calc_kappa.MoyoDataset") as m_moyo,
        patch("matbench_discovery.phonons.calc_kappa.MoyoAdapter"),
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
    from matbench_discovery.phonons.calc_kappa import calc_kappa_for_structure

    freqs = (
        np.array([[-1.0, 0.1, 0.2, 1.0, 2.0, 3.0]])
        if has_imag
        else np.array([[0.001, 0.002, 0.003, 1.0, 2.0, 3.0]])
    )

    with (
        patch("matbench_discovery.phonons.calc_kappa.MoyoDataset") as m_moyo,
        patch("matbench_discovery.phonons.calc_kappa.MoyoAdapter"),
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
    from matbench_discovery.phonons.calc_kappa import calc_kappa_for_structure

    freqs_2d = np.array([[0.001, 0.002, 0.003, 1.0, 2.0, 3.0]])
    with (
        patch("matbench_discovery.phonons.calc_kappa.MoyoDataset") as m_moyo,
        patch("matbench_discovery.phonons.calc_kappa.MoyoAdapter"),
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
    from matbench_discovery.phonons.calc_kappa import calc_kappa_for_structure

    freqs_2d = np.array([[0.001, 0.002, 0.003, 1.0, 2.0, 3.0]])
    with (
        patch("matbench_discovery.phonons.calc_kappa.MoyoDataset") as m_moyo,
        patch("matbench_discovery.phonons.calc_kappa.MoyoAdapter"),
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
