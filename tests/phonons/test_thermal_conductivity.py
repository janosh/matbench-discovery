"""Tests for thermal conductivity calculation module."""

import itertools

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk
from ase.calculators.calculator import Calculator
from ase.calculators.emt import EMT
from phono3py.api_phono3py import Phono3py
from phono3py.conductivity.calculators import RTACalculator
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from pymatviz.enums import Key

from matbench_discovery.enums import MbdKey
from matbench_discovery.phonons import thermal_conductivity as ltc
from matbench_discovery.phonons.thermal_conductivity import calculate_fc2_set


@pytest.mark.parametrize(
    "forces",
    [
        np.zeros((2, 2)),
        np.full((2, 3), np.nan),
        np.full((2, 3), 1j),
    ],
)
def test_force_validation_rejects_malformed_predictions(forces: np.ndarray) -> None:
    """Force predictions require exact shape, real dtype, and finite values."""
    with pytest.raises(ValueError, match="Invalid test forces"):
        ltc._validate_forces(forces, (2, 3), "test")  # noqa: SLF001


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

    custom_mesh = (4, 4, 4)
    ph3 = ltc.init_phono3py(
        test_atoms,
        fc2_supercell=fc2_supercell,
        fc3_supercell=fc3_supercell,
        q_point_mesh=custom_mesh,
    )
    assert ph3.mesh_numbers is not None
    assert tuple(ph3.mesh_numbers.tolist()) == custom_mesh

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


def test_calculate_fc3_set(test_ph3: Phono3py, test_calculator: EMT) -> None:
    """Test calculation of 3rd order force constants."""
    forces = ltc.calculate_fc3_set(test_ph3, test_calculator, pbar_kwargs={})
    assert isinstance(forces, np.ndarray)
    assert forces.ndim == 3
    assert forces.shape[-1] == 3

    capped_forces = ltc.calculate_fc3_set(
        test_ph3,
        test_calculator,
        pbar_kwargs={},
        max_evaluations=0,
    )
    assert capped_forces.shape == forces.shape
    assert np.count_nonzero(capped_forces) == 0


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


def test_calculate_conductivity(test_ph3: Phono3py, test_calculator: EMT) -> None:
    """Test thermal conductivity calculation."""
    # First need to compute force constants
    test_ph3, fc2_set, _ = ltc.get_fc2_and_freqs(
        test_ph3, test_calculator, pbar_kwargs={}
    )

    fc3_set = ltc.calculate_fc3_set(test_ph3, calculator=test_calculator)

    test_ph3.phonon_forces = fc2_set
    test_ph3.forces = fc3_set
    test_ph3.produce_fc2(symmetrize_fc2=True)
    test_ph3.produce_fc3(symmetrize_fc3r=True)

    ph3, kappa_dict, kappa = ltc.calculate_conductivity(test_ph3, [300])

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
    """Coherence term reduces to a plain band sum for constant heat capacity."""
    rng = np.random.default_rng(seed=0)
    # array axes are T, q-points, bands, xyz (coherence has an extra bands axis)
    mode_kappa_p_rta = rng.random((2, 3, 3, 3))
    mode_kappa_coherence = rng.random((2, 3, 3, 3, 3))
    heat_capacity = np.full((2, 3, 3), 0.7)

    # with C_i == C_j the weight 2 * C_i / (C_i + C_j) is exactly 1
    result = ltc.calc_mode_kappa_tot(
        mode_kappa_p_rta, mode_kappa_coherence, heat_capacity
    )
    np.testing.assert_allclose(
        result, mode_kappa_p_rta + mode_kappa_coherence.sum(axis=2), rtol=1e-12, atol=0
    )

    # zero heat capacity yields 0/0 = NaN which is zeroed, leaving only the RTA term
    result_zero_cap = ltc.calc_mode_kappa_tot(
        mode_kappa_p_rta, mode_kappa_coherence, np.zeros_like(heat_capacity)
    )
    np.testing.assert_array_equal(result_zero_cap, mode_kappa_p_rta)


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
    """FC2 evaluation preserves the calculator forces and displacement shape."""
    atoms = make_si2_phonopy_atoms()
    ph3 = build_ph3_with_fc2(atoms, np.eye(3, dtype=int), distance=0.03)
    expected_forces = np.hstack((0.1 * np.eye(3), -0.1 * np.eye(3)))[0].reshape(2, 3)
    calc = MockCalculator(expected_forces)
    force_set = calculate_fc2_set(ph3, calc, pbar_kwargs={"disable": True})
    assert force_set.shape == (
        len(ph3.phonon_supercells_with_displacements),
        len(ph3.phonon_supercell),
        3,
    )
    np.testing.assert_allclose(
        force_set, np.broadcast_to(expected_forces, force_set.shape)
    )


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


def test_harmonic_data_from_phono3py(test_ph3: Phono3py, test_calculator: EMT) -> None:
    """FC2 yields a seekpath band structure with normalized eigenvectors, DOS and Cv."""
    from matbench_discovery.phonons import harmonic

    ph3, _fc2_set, mesh_freqs = ltc.get_fc2_and_freqs(
        test_ph3, test_calculator, pbar_kwargs={"disable": True}
    )
    data = harmonic.harmonic_data_from_phono3py(ph3)
    n_atoms = len(ph3.primitive)  # Cu FCC primitive has 1 atom
    assert n_atoms == 1
    band_path = data["band_path"]
    n_q = len(band_path["q_points"])
    segments = band_path["segments"]
    # seekpath's FCC path is GAMMA-X-U|K-GAMMA-L-W-X: connected segments share their
    # endpoint q-point while the U|K break advances by one
    labels = [(seg["start_label"], seg["end_label"]) for seg in segments]
    assert labels == [
        ("GAMMA", "X"),
        ("X", "U"),
        ("K", "GAMMA"),
        ("GAMMA", "L"),
        ("L", "W"),
        ("W", "X"),
    ]
    assert segments[0]["start_index"] == 0
    assert segments[-1]["end_index"] == n_q - 1
    for previous, current in itertools.pairwise(segments):
        shared = previous["end_label"] == current["start_label"]
        assert current["start_index"] == previous["end_index"] + (not shared)
    # q-points are spaced by (just over) the reference distance along each segment
    distances = band_path["distances"]
    start, end = segments[0]["start_index"], segments[0]["end_index"]
    spacing = (distances[end] - distances[start]) / (end - start)
    assert 1 <= spacing / harmonic.BAND_PATH_REFERENCE_DISTANCE <= 1.05
    assert band_path["q_points"].shape == (n_q, 3)
    assert band_path["distances"].shape == (n_q,)
    assert band_path["frequencies"].shape == (n_q, 3 * n_atoms)
    assert band_path["group_velocities"].shape == (n_q, 3 * n_atoms, 3)
    eigenvectors = band_path["eigenvectors"]
    assert eigenvectors.shape == (n_q, 3 * n_atoms, n_atoms, 3, 2)
    # stored to EIGENVECTOR_DECIMALS: rounding each of the n = 6 * n_atoms components
    # by at most 0.5 * 10^-d shifts |e|^2 by at most 2 * 0.5 * 10^-d * sum|e_i|
    # <= sqrt(n) * 10^-d (Cauchy-Schwarz with |e| = 1)
    np.testing.assert_array_equal(
        eigenvectors, np.round(eigenvectors, harmonic.EIGENVECTOR_DECIMALS)
    )
    np.testing.assert_allclose(
        eigenvectors.reshape(n_q, 3 * n_atoms, -1).__pow__(2).sum(axis=-1),
        1,
        rtol=0,
        atol=np.sqrt(6 * n_atoms) * 10**-harmonic.EIGENVECTOR_DECIMALS,
    )
    # Gamma point: acoustic modes at zero frequency, zero distance
    np.testing.assert_allclose(band_path["q_points"][0], 0, atol=0)
    assert band_path["distances"][0] == 0
    np.testing.assert_allclose(band_path["frequencies"][0], 0, rtol=0, atol=1e-6)
    # band frequencies live in the same window as the phono3py mesh frequencies
    assert band_path["frequencies"].max() == pytest.approx(mesh_freqs.max(), rel=0.05)
    dos = data["dos"]
    assert dos["mesh"] == [2, 2, 2]
    assert len(dos["frequencies"]) == len(dos["densities"])
    # the coarse fixture mesh under-resolves the tetrahedron DOS (2.25 states); on a
    # converged mesh the DOS integrates to 3N states per primitive cell
    dense_dos = harmonic.harmonic_phonon_data(
        harmonic.phonopy_from_phono3py(ph3), mesh=(12, 12, 12)
    )["dos"]
    n_states = np.trapezoid(dense_dos["densities"], dense_dos["frequencies"])
    assert n_states == pytest.approx(3 * n_atoms, rel=1e-3)
    thermal = data["thermal_properties"]
    assert thermal["temperatures"][0] == 0
    assert thermal["temperatures"][-1] == harmonic.THERMAL_T_MAX
    assert thermal["heat_capacity"][0] == 0
    # high-T heat capacity approaches the classical 3N k_B = 24.94 J/K/mol per atom
    assert thermal["heat_capacity"][-1] == pytest.approx(3 * n_atoms * 8.314, rel=0.02)
    primitive = data["primitive"]
    assert primitive["symbols"] == ["Cu"]
    assert primitive["lattice"].shape == (3, 3)
    assert primitive["frac_coords"].shape == (n_atoms, 3)
    assert data["frequency_unit"] == "THz"


def test_seekpath_band_path_transfers_q_points_to_phonopy_basis(
    test_atoms: Atoms, test_ph3: Phono3py, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Path q-points land on the intended Cartesian k even when seekpath re-bases."""
    import seekpath

    from matbench_discovery.phonons import harmonic

    def emt_phonon(atoms: Atoms) -> Phonopy:
        ph3 = ltc.init_phono3py(
            atoms,
            fc2_supercell=test_atoms.info["fc2_supercell"],
            fc3_supercell=test_atoms.info["fc3_supercell"],
            q_point_mesh=test_atoms.info["q_point_mesh"],
        )
        return harmonic.phonopy_from_phono3py(
            ltc.get_fc2_and_freqs(ph3, EMT(), pbar_kwargs={"disable": True})[0]
        )

    pristine = harmonic.phonopy_from_phono3py(
        ltc.get_fc2_and_freqs(test_ph3, EMT(), pbar_kwargs={"disable": True})[0]
    )
    assert len(harmonic.seekpath_band_path(pristine)["explicit_segments"]) == 6

    # a 0.3 % xy shear breaks FCC down to C2/m: seekpath's standardized primitive is
    # then a different basis from phonopy's, so the naive fractional transfer would
    # sample the wrong k-points (several THz off), while the Cartesian round trip
    # stays within the strain-induced shift of the pristine cell's frequencies
    sheared_atoms = test_atoms.copy()
    shear = np.eye(3)
    shear[0, 1] = 0.003
    sheared_atoms.set_cell(np.asarray(test_atoms.cell) @ shear, scale_atoms=True)
    sheared = emt_phonon(sheared_atoms)
    path = harmonic.seekpath_band_path(sheared)
    assert path["spacegroup_international"] == "C2/m"
    sheared.run_qpoints(path["explicit_kpoints_rel"])
    sheared_freqs = sheared.get_qpoints_dict()["frequencies"]
    k_cart = path["explicit_kpoints_rel"] @ np.linalg.inv(sheared.primitive.cell).T
    pristine.run_qpoints(k_cart @ np.asarray(pristine.primitive.cell).T)
    pristine_freqs = pristine.get_qpoints_dict()["frequencies"]
    np.testing.assert_allclose(sheared_freqs, pristine_freqs, rtol=0, atol=0.1)
    naive_q = seekpath.get_explicit_k_path(
        sheared.primitive.totuple(),
        reference_distance=harmonic.BAND_PATH_REFERENCE_DISTANCE,
    )["explicit_kpoints_rel"]
    sheared.run_qpoints(naive_q)
    assert np.abs(sheared.get_qpoints_dict()["frequencies"] - pristine_freqs).max() > 1

    # a seekpath primitive that is not the same lattice at all still aborts
    real_get_path = seekpath.get_explicit_k_path

    def stretched_get_path(*args: object, **kwargs: object) -> dict[str, object]:
        path = dict(real_get_path(*args, **kwargs))
        path["primitive_lattice"] = np.asarray(path["primitive_lattice"]) * 1.01
        return path

    monkeypatch.setattr(seekpath, "get_explicit_k_path", stretched_get_path)
    with pytest.raises(ValueError, match="not a re-based phonopy primitive"):
        harmonic.seekpath_band_path(pristine)


def test_harmonic_eigenvectors_are_dynamical_matrix_eigenvectors() -> None:
    """Stored [band][atom][xyz] eigenvectors solve D e = w^2 e for a 2-atom cell."""
    from matbench_discovery.phonons import harmonic

    atoms = bulk("CuAu", "cesiumchloride", a=3.0)
    ph3 = ltc.init_phono3py(
        atoms,
        fc2_supercell=2 * np.eye(3),
        fc3_supercell=np.eye(3),
        q_point_mesh=(2, 2, 2),
    )
    phonon = harmonic.phonopy_from_phono3py(
        ltc.get_fc2_and_freqs(ph3, EMT(), pbar_kwargs={"disable": True})[0]
    )
    band_path = harmonic.harmonic_phonon_data(phonon, mesh=(2, 2, 2))["band_path"]
    n_atoms = len(phonon.primitive)
    assert n_atoms == 2
    # a transposed (wrong) layout would still pass the normalization check because
    # the eigenvector matrix is unitary, so verify the eigen-equation directly
    for q_idx in range(0, len(band_path["q_points"]), 20):
        phonon.dynamical_matrix.run(band_path["q_points"][q_idx])
        dyn_mat = phonon.dynamical_matrix.dynamical_matrix
        stored = band_path["eigenvectors"][q_idx]
        eig_vecs = (stored[..., 0] + 1j * stored[..., 1]).reshape(3 * n_atoms, -1)
        freqs = band_path["frequencies"][q_idx]
        eig_vals = np.sign(freqs) * (freqs / phonon.unit_conversion_factor) ** 2
        for band_idx, eig_vec in enumerate(eig_vecs):
            np.testing.assert_allclose(
                dyn_mat @ eig_vec, eig_vals[band_idx] * eig_vec, rtol=0, atol=1e-4
            )
    # at Gamma the acoustic displacement is uniform: |e_atom| scales with sqrt(mass)
    gamma_acoustic = band_path["eigenvectors"][0, :3]
    per_atom_norm = np.sqrt((gamma_acoustic**2).sum(axis=(2, 3)))  # (3 bands, 2 atoms)
    displacement = per_atom_norm / np.sqrt(np.asarray(phonon.primitive.masses))
    np.testing.assert_allclose(displacement[:, 1], displacement[:, 0], rtol=1e-3)
