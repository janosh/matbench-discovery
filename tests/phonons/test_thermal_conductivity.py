"""Tests for force-constant, conductivity and harmonic-phonon helpers."""

import itertools
from collections.abc import Sequence

import numpy as np
import pytest
from matscikit import EMT, Structure
from pymatviz.enums import Key

from matbench_discovery.enums import MbdKey
from matbench_discovery.phonons import harmonic
from matbench_discovery.phonons import thermal_conductivity as ltc
from matbench_discovery.phonons.adapters import StandardKappaAdapter
from matbench_discovery.phonons.pipeline import KappaSettings

ADAPTER = StandardKappaAdapter()
SETTINGS = KappaSettings()


@pytest.fixture
def cu_atoms() -> Structure:
    """FCC Cu primitive cell (supported by EMT) with PhononDB-style metadata."""
    a = 3.61
    lattice = a / 2 * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]])
    return Structure(
        lattice,
        ["Cu"],
        [[0, 0, 0]],
        fractional=True,
        # a 2x2x2 fc3 supercell of this 1-atom cell makes every atom its own inversion
        # image about atom 0, forcing fc3 = 0 by symmetry; 3x3x3 keeps it finite
        info={
            "fc2_supercell": np.diag([2, 2, 2]),
            "fc3_supercell": np.diag([3, 3, 3]),
            "q_point_mesh": (4, 4, 4),
        },
    )


@pytest.fixture
def cu_calc(cu_atoms: Structure) -> ltc.KappaCalculation:
    """Phonons with displacements on the 2x2x2 benchmark-style setup."""
    return ltc.init_phonons(
        cu_atoms,
        fc2_supercell=cu_atoms.info["fc2_supercell"],
        fc3_supercell=cu_atoms.info["fc3_supercell"],
        q_point_mesh=cu_atoms.info["q_point_mesh"],
        displacement_distance=0.01,  # golden values below were computed at 0.01 Å
    )


def test_init_phonons(cu_calc: ltc.KappaCalculation) -> None:
    """Setup builds both supercells and symmetry-reduced displacement sets."""
    assert cu_calc.mesh == (4, 4, 4)
    assert len(cu_calc.phonons.supercell) == 27
    assert len(cu_calc.phonons.phonon_supercell) == 8
    assert len(cu_calc.phonons.phonon_displacements) == 1  # one fcc displacement
    configs = cu_calc.phonons.displacements
    assert len(configs[0]) == 1  # singles first, then pairs
    assert all(len(config) in (1, 2) for config in configs)
    # default step matches the 0.03 Å displacements of the PhononDB-PBE reference
    default_calc = ltc.init_phonons(
        cu_calc.phonons.unitcell, fc2_supercell=np.eye(3), fc3_supercell=np.eye(3)
    )
    _, fc2_displacement = default_calc.phonons.phonon_displacements[0]
    assert np.linalg.norm(fc2_displacement) == pytest.approx(0.03, abs=1e-12)
    whole_floats = ltc.init_phonons(
        cu_calc.phonons.unitcell,
        fc2_supercell=np.eye(3),
        fc3_supercell=np.eye(3),
        q_point_mesh=(4.0, 4, 4),  # ty: ignore[invalid-argument-type]
    )
    assert whole_floats.mesh == (4, 4, 4)
    for bad_mesh in [(2, 2), (4.9, 4, 4), (0, 4, 4), (np.nan, 4, 4)]:
        with pytest.raises(ValueError, match="3 positive integers"):
            ltc.init_phonons(
                cu_calc.phonons.unitcell,
                fc2_supercell=np.eye(3),
                fc3_supercell=np.eye(3),
                q_point_mesh=bad_mesh,  # ty: ignore[invalid-argument-type]
            )


class CountingModel:
    """EMT wrapper recording the size of every model call."""

    def __init__(self, forces_scale: float = 1.0) -> None:
        """Start with no recorded calls; ``forces_scale`` multiplies EMT's forces."""
        self.batch_sizes: list[int] = []
        self.forces_scale = forces_scale

    def compute(
        self, structures: Sequence[Structure], properties: Sequence[str]
    ) -> list[dict]:
        """Delegate to EMT after recording the batch size."""
        self.batch_sizes.append(len(structures))
        results = EMT().compute(structures, properties)
        return [{**res, "forces": res["forces"] * self.forces_scale} for res in results]


def test_fc3_forces_batch_and_cap(cu_calc: ltc.KappaCalculation) -> None:
    """Forces are evaluated in batches; capped evaluations leave zero slots."""
    model = CountingModel()
    ADAPTER.calculate_fc3(cu_calc, model, KappaSettings(batch_size=3))
    forces = cu_calc.phonons.fc3_forces
    assert forces is not None
    assert forces.shape == (len(cu_calc.phonons.displacements), 27, 3)
    assert model.batch_sizes[:-1] == [3] * (len(model.batch_sizes) - 1)
    assert sum(model.batch_sizes) == len(forces)
    # matscikit's native EMT gives the same forces bit for bit, batched or not
    ADAPTER.calculate_fc3(cu_calc, EMT(), SETTINGS)
    np.testing.assert_array_equal(cu_calc.phonons.fc3_forces, forces)
    ADAPTER.calculate_fc3(cu_calc, EMT(), SETTINGS, max_evaluations=2)
    capped = cu_calc.phonons.fc3_forces
    assert capped is not None
    np.testing.assert_array_equal(capped[:2], forces[:2])
    assert not capped[2:].any()
    with pytest.raises(ValueError, match="supercell 0 are not finite"):
        ADAPTER.calculate_fc3(cu_calc, CountingModel(forces_scale=np.nan), SETTINGS)


def test_calculate_conductivity(cu_calc: ltc.KappaCalculation) -> None:
    """End-to-end EMT Cu kappa: finite, isotropic and consistent with mode sums."""
    freqs = ADAPTER.calculate_fc2(cu_calc, EMT(), SETTINGS)
    assert cu_calc.phonons.fc2_forces is not None
    assert cu_calc.phonons.fc2_forces.shape == (1, 8, 3)
    assert freqs.shape == (64, 3)  # every point of the 4x4x4 grid, 3 bands
    np.testing.assert_allclose(freqs[0], 0, atol=1e-5)  # acoustic modes at Gamma
    # grid index a0 + n0 (a1 + n1 a2): equal to diagonalizing every grid point directly
    # (star members differ at eigensolver round-off, measured 3.6e-15 THz)
    address = np.indices(cu_calc.mesh[::-1]).reshape(3, -1).T[:, ::-1]
    direct = cu_calc.phonons.harmonic.frequencies(address / np.asarray(cu_calc.mesh))
    np.testing.assert_allclose(freqs, direct, rtol=0, atol=1e-12)
    # golden values pin matscikit numerics (fc2 fit, dynamical matrix, fc3, Wigner RTA);
    # EMT and the algorithms are deterministic, so 1e-9 leaves only summation-order
    # round-off: an intentional numerics change must update these
    np.testing.assert_allclose(
        freqs[1], [2.42804731, 2.42804731, 5.45674688], rtol=1e-8
    )
    ADAPTER.calculate_fc3(cu_calc, EMT(), SETTINGS)
    assert cu_calc.phonons.fc3_forces is not None
    assert cu_calc.phonons.fc3_forces.shape == (
        len(cu_calc.phonons.displacements),
        27,
        3,
    )
    kappa = ltc.calculate_conductivity(cu_calc, [300])
    kappa_tot = kappa[MbdKey.kappa_tot_rta]
    assert kappa_tot.shape == (1, 6)
    # matscikit 61a0f69 (bit-identical on 1, 2, 4 and 14 rayon threads)
    np.testing.assert_allclose(kappa_tot[0, :3], 22.107493045885906, rtol=1e-9)
    np.testing.assert_allclose(kappa_tot[0, 3:], 0, atol=1e-12)
    np.testing.assert_allclose(
        kappa_tot, kappa[MbdKey.kappa_p_rta] + kappa[MbdKey.kappa_c], rtol=1e-12
    )
    weights = kappa[Key.mode_weights]
    assert weights.sum() == 64
    assert kappa[Key.q_points].shape == (len(weights), 3)
    assert kappa[Key.ph_freqs].shape == (len(weights), 3)
    mode_kappa = kappa[MbdKey.mode_kappa_tot_rta]
    assert mode_kappa.shape == (1, len(weights), 3, 6)
    # matscikit sums its per-mode split exactly (measured relative gap 3.2e-16)
    np.testing.assert_allclose(
        mode_kappa.sum(axis=(1, 2)) / weights.sum(), kappa_tot, rtol=1e-13, atol=1e-14
    )


def test_harmonic_data(cu_calc: ltc.KappaCalculation) -> None:
    """FC2 yields an HPKOT band structure with normalized eigenvectors, DOS and Cv."""
    mesh_freqs = ADAPTER.calculate_fc2(cu_calc, EMT(), SETTINGS)
    data = ADAPTER.harmonic_phonons(cu_calc, SETTINGS)
    assert data["errors"] == {}
    mesh = data["mesh"]
    assert np.asarray(mesh["weights"]).sum() == np.prod(mesh["numbers"])
    assert mesh["frequencies"].shape[0] == len(mesh["q_points"]) == len(mesh["weights"])
    # stable Cu: acoustic modes at Gamma must not be reported as imaginary
    assert data["thermal_properties"]["n_imaginary_modes_excluded"] == 0
    n_atoms = len(cu_calc.phonons.primitive)  # Cu FCC primitive has 1 atom
    assert n_atoms == 1
    band_path = data["band_path"]
    n_q = len(band_path["q_points"])
    segments = band_path["segments"]
    # HPKOT's cF2 path is GAMMA-X-U|K-GAMMA-L-W-X: connected segments share their
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
    assert band_path["frequencies"].shape == (n_q, 3 * n_atoms)
    assert band_path["group_velocities"].shape == (n_q, 3 * n_atoms, 3)
    eigenvectors = band_path["eigenvectors"]
    assert eigenvectors.shape == (n_q, 3 * n_atoms, n_atoms, 3, 2)
    # stored to EIGENVECTOR_DECIMALS: rounding each of the n = 6 * n_atoms components
    # by at most 0.5 * 10^-d shifts |e|^2 by at most sqrt(n) * 10^-d
    np.testing.assert_allclose(
        eigenvectors.reshape(n_q, 3 * n_atoms, -1).__pow__(2).sum(axis=-1),
        1,
        rtol=0,
        atol=np.sqrt(6 * n_atoms) * 10**-harmonic.EIGENVECTOR_DECIMALS,
    )
    # Gamma point: acoustic modes at zero frequency, zero distance
    np.testing.assert_allclose(band_path["q_points"][0], 0, atol=0)
    assert band_path["distances"][0] == 0
    np.testing.assert_allclose(band_path["frequencies"][0], 0, rtol=0, atol=1e-5)
    assert band_path["frequencies"].max() == pytest.approx(mesh_freqs.max(), rel=0.05)
    dos = data["dos"]
    assert dos["mesh"] == [4, 4, 4]
    assert len(dos["frequencies"]) == len(dos["densities"])
    # on a converged mesh the tetrahedron DOS integrates to 3N states per primitive cell
    dense_dos = harmonic.harmonic_phonon_data(
        cu_calc.phonons.harmonic, mesh=(12, 12, 12)
    )["dos"]
    n_states = np.trapezoid(dense_dos["densities"], dense_dos["frequencies"])
    assert n_states == pytest.approx(3 * n_atoms, rel=1e-3)
    thermal = data["thermal_properties"]
    assert thermal["temperatures"][0] == 0
    assert thermal["temperatures"][-1] == harmonic.THERMAL_TEMPERATURES[-1] == 1000
    assert thermal["heat_capacity"][0] == 0
    # high-T heat capacity approaches 3N k_B = 24.94 J/K/mol per atom times the
    # fraction of included modes: the 3 acoustic modes at Gamma (of 64 q x 3 bands on
    # the 4x4x4 mesh) are excluded from thermodynamic sums
    included = 1 - 3 / (np.prod(mesh["numbers"]) * 3 * n_atoms)
    expected_cv = 3 * n_atoms * 8.314 * included
    assert thermal["heat_capacity"][-1] == pytest.approx(expected_cv, rel=0.02)
    primitive = data["primitive"]
    assert primitive["symbols"] == ["Cu"]
    assert primitive["lattice"].shape == (3, 3)
    assert data["frequency_unit"] == "THz"


def test_harmonic_artifact_round_trip(cu_calc: ltc.KappaCalculation) -> None:
    """Stored FC2 restores the same phonons; a foreign schema is rejected."""
    ADAPTER.calculate_fc2(cu_calc, EMT(), SETTINGS)
    data = harmonic.harmonic_phonon_data(cu_calc.phonons.harmonic, mesh=(2, 2, 2))
    restored = harmonic.phonons_from_harmonic_data(data)
    # compact FC2 is re-expanded by the lattice translations; the symmetrized fit is
    # translation invariant to a few ulp of the largest force constant (2.2e-16)
    fc2 = cu_calc.phonons.harmonic.force_constants
    fc2_atol = 8 * np.finfo(np.float64).eps * np.abs(fc2).max()
    np.testing.assert_allclose(restored.force_constants, fc2, rtol=0, atol=fc2_atol)
    q_points = np.array([[0.1, 0.2, 0.3], [0.5, 0, 0]])
    np.testing.assert_allclose(
        restored.frequencies(q_points),
        cu_calc.phonons.harmonic.frequencies(q_points),
        rtol=1e-12,
        atol=1e-12,
    )
    with pytest.raises(ValueError, match="Unsupported harmonic schema"):
        harmonic.phonons_from_harmonic_data({**data, "schema_version": -1})
    shuffled = {**data, "fc2": {**data["fc2"]}}
    shuffled["fc2"]["supercell_frac_coords"] = data["fc2"]["supercell_frac_coords"][
        ::-1
    ]
    with pytest.raises(ValueError, match="supercell order"):
        harmonic.phonons_from_harmonic_data(shuffled)


class FailingPhonons:
    """Delegate to real Phonons but raise in one analysis method."""

    def __init__(self, phonons: object, method: str) -> None:
        """Wrap ``phonons``; calls to ``method`` raise RuntimeError."""
        self._phonons, self._method = phonons, method

    def __getattr__(self, name: str) -> object:
        """Forward attributes, failing the configured method."""
        if name == self._method:

            def fail(*_args: object, **_kwargs: object) -> None:
                raise RuntimeError(f"failed {name}")

            return fail
        return getattr(self._phonons, name)


@pytest.mark.parametrize(
    ("method", "section", "absent"),
    [
        ("modes", "band_path", {"band_path"}),
        # DOS and thermal properties come from the mesh's diagonalization
        ("mesh", "mesh", {"mesh", "dos", "thermal_properties"}),
    ],
)
def test_harmonic_analysis_failures_are_independent(
    cu_calc: ltc.KappaCalculation, method: str, section: str, absent: set[str]
) -> None:
    """An analysis failure keeps FC2 and every independent successful analysis."""
    ADAPTER.calculate_fc2(cu_calc, EMT(), SETTINGS)
    data = harmonic.harmonic_phonon_data(
        FailingPhonons(cu_calc.phonons.harmonic, method),  # ty: ignore[invalid-argument-type]
        mesh=(2, 2, 2),
    )
    assert set(data["errors"]) == {section}
    assert data["errors"][section]["message"] == f"RuntimeError: failed {method}"
    assert "Traceback" in data["errors"][section]["traceback"]
    assert "fc2" in data
    sections = {"band_path", "mesh", "dos", "thermal_properties"}
    assert sections & data.keys() == sections - absent


def test_harmonic_eigenvectors_solve_gamma_dynamical_matrix() -> None:
    """Stored [band][atom][xyz] eigenvectors solve an independent D(Γ) e = w² e."""
    from matscikit import Phonons

    cell = Structure(
        3.0 * np.eye(3), ["Cu", "Au"], [[0, 0, 0], [0.5, 0.5, 0.5]], fractional=True
    )
    phonons = Phonons(cell, [2, 2, 2])
    phonons.generate_displacements(distance=0.01)
    forces = [
        result["forces"]
        for result in EMT().compute(phonons.supercells_with_displacements(), ["forces"])
    ]
    phonons.produce_fc2(np.array(forces))
    band_path = harmonic.harmonic_band_data(phonons)
    assert np.allclose(band_path["q_points"][0], 0)  # HPKOT paths start at Gamma
    # D(Γ)_{pa,qb} = Σ_{j -> q} Φ(rep(p), j)_{ab} / sqrt(m_p m_q): no phases at Γ, so a
    # transposed or misordered eigenvector layout cannot satisfy it
    supercell, primitive = phonons.supercell, phonons.primitive
    frac_in_prim = supercell.positions @ np.linalg.inv(primitive.lattice)
    diff = frac_in_prim[:, None] - primitive.frac_coords[None]
    s2p = np.argmin(np.linalg.norm(diff - np.round(diff), axis=-1), axis=1)
    reps = [int(np.flatnonzero(s2p == p)[0]) for p in range(len(primitive))]
    fc2, masses = phonons.force_constants, primitive.masses
    n_prim = len(primitive)
    dyn = np.zeros((3 * n_prim, 3 * n_prim))
    for p, rep in enumerate(reps):
        for j in range(len(supercell)):
            q = s2p[j]
            dyn[3 * p : 3 * p + 3, 3 * q : 3 * q + 3] += fc2[rep, j] / np.sqrt(
                masses[p] * masses[q]
            )
    stored = band_path["eigenvectors"][0]
    eig_vecs = (stored[..., 0] + 1j * stored[..., 1]).reshape(3 * n_prim, 3 * n_prim)
    freqs = band_path["frequencies"][0]
    # THz per sqrt(eV/A^2/amu) with matscikit's CODATA 2014 constants
    factor = np.sqrt(1.6021766208e-19 / (1e-20 * 1.660539040e-27)) / (2 * np.pi * 1e12)
    eig_vals = np.sign(freqs) * (freqs / factor) ** 2
    for band_idx, eig_vec in enumerate(eig_vecs):
        np.testing.assert_allclose(
            dyn @ eig_vec, eig_vals[band_idx] * eig_vec, rtol=0, atol=1e-4
        )


def test_harmonic_thermal_properties_report_excluded_imaginary_modes() -> None:
    """Thermal properties count modes excluded below the imaginary threshold."""
    from matscikit import Phonons

    cell = Structure(
        4.0 * np.eye(3), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]], fractional=True
    )
    phonons = Phonons(cell, [1, 1, 1])
    rng = np.random.default_rng(seed=0)
    n_atoms = len(phonons.supercell)
    force_constants = rng.normal(scale=0.5, size=(n_atoms, n_atoms, 3, 3))
    # symmetric FC2 -> real eigenvalues, a mix of stable and imaginary branches
    phonons.set_force_constants(force_constants + force_constants.transpose(1, 0, 3, 2))
    thermal = harmonic.harmonic_thermal_data(
        phonons.mesh([4, 4, 4], temperatures_k=harmonic.THERMAL_TEMPERATURES.tolist())
    )
    assert 0 < thermal["n_imaginary_modes_excluded"] <= thermal["n_mesh_modes"]
    # the dropped modes never contribute: Cv stays below the classical 3 N k_B limit
    assert np.asarray(thermal["heat_capacity"])[-1] < 3 * len(phonons.primitive) * 8.314
