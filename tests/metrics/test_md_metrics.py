"""Tests for molecular dynamics trajectory metrics."""

import itertools
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from ase import Atoms, units
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.singlepoint import SinglePointCalculator
from ase.data import atomic_masses

from matbench_discovery.enums import Model
from matbench_discovery.metrics import md as md_metrics

np_rng = np.random.default_rng(seed=0)
ANISOTROPIC_STRESS_VOIGT = np.array(
    [-1.0 * units.GPa, -2.0 * units.GPa, -6.0 * units.GPa, 0, 0, 0]
)


class ConstantCalculator(Calculator):
    """Predicts a fixed energy and uniform force regardless of geometry. Geometry-
    independent so it works with the array-based pipeline (which feeds the calculator
    fresh Atoms built from positions only, without per-frame info/arrays).
    """

    implemented_properties = ("energy", "forces")

    def __init__(self, energy: float = 0, force: float = 0, **kwargs: object) -> None:
        """Store the constant energy and per-component force to return."""
        super().__init__(**kwargs)
        self.energy = energy
        self.force = force

    def calculate(
        self,
        atoms: Atoms | None = None,
        properties: list[str] | None = None,
        system_changes: list[str] = all_changes,
    ) -> None:
        """Return the stored constant energy and a uniform force on every atom."""
        super().calculate(atoms, properties, system_changes)
        assert atoms is not None
        self.results["energy"] = self.energy
        self.results["forces"] = np.full((len(atoms), 3), self.force)


class BrokenStressCalculator(Calculator):
    """Calculator that raises an unexpected RuntimeError for stress."""

    implemented_properties = ("stress",)

    def calculate(
        self,
        atoms: Atoms | None = None,
        properties: list[str] | None = None,
        system_changes: list[str] = all_changes,
    ) -> None:
        """Raise a backend-like failure when stress is requested."""
        super().calculate(atoms, properties, system_changes)
        raise RuntimeError("simulated stress backend failure")


def attach_ref_labels(atoms: Atoms, energy: float, forces: np.ndarray) -> Atoms:
    """Attach reference energy/forces as SinglePointCalculator results so
    Trajectory.from_ase extracts them as the reference labels.
    """
    atoms.calc = SinglePointCalculator(atoms, energy=energy, forces=forces)
    return atoms


def make_jiggled_frames(
    n_frames: int, *, pressure_gpa: float | None = None
) -> list[Atoms]:
    """Periodic 8-atom cubic-lattice frames with small random displacements. The 2.4 A
    spacing keeps Si-Si pairs within the covalent-radius bond cutoff so calc_adf finds
    bonded triplets (real Si-Si bond ~2.35 A; covalent cutoff ~2.78 A).
    """
    base_positions = 2.4 * np.array(
        list(itertools.product(range(2), repeat=3)), dtype=float
    )
    frames = []
    for _ in range(n_frames):
        positions = base_positions + np_rng.normal(scale=0.05, size=(8, 3))
        atoms = Atoms("Si8", positions=positions, cell=[6] * 3, pbc=True)
        if pressure_gpa is not None:
            stress = -pressure_gpa * units.GPa * np.eye(3)
            atoms.calc = SinglePointCalculator(atoms, stress=stress)
        frames.append(atoms)
    return frames


def h2_frame(
    energy: float,
    *,
    positions: np.ndarray | None = None,
    forces: np.ndarray | None = None,
) -> Atoms:
    """An H2 reference frame (random positions unless given) carrying energy/forces."""
    if positions is None:
        positions = np_rng.random((2, 3)) * 2
    if forces is None:
        forces = np.zeros((2, 3))
    atoms = Atoms("H2", positions=positions, cell=[5] * 3)
    return attach_ref_labels(atoms, energy=energy, forces=forces)


# === RDF ===


# fixed cell, and alternating very different volumes (probes per-frame densities)
@pytest.mark.parametrize("box_lens", [(10,), (10, 16)])
def test_calc_rdf_ideal_gas_is_flat(box_lens: tuple[float, ...]) -> None:
    """g(r) of uniformly random positions should match paper Eq. 4 normalization."""
    rng = np.random.default_rng(seed=0)
    n_atoms = 64
    frames = [
        Atoms(
            "Ar64", positions=rng.random((n_atoms, 3)) * box, cell=[box] * 3, pbc=True
        )
        for _ in range(40)
        for box in box_lens
    ]
    radii, g_r = md_metrics.calc_rdf(frames, n_bins=50)
    # Paper Eq. 4 uses rho*N, so finite-N ideal gas plateaus at (N - 1) / N.
    np.testing.assert_allclose(
        g_r[radii > 2].mean(), (n_atoms - 1) / n_atoms, atol=0.01
    )


def test_calc_rdf_peak_position() -> None:
    """g(r) of a near-rigid dimer should peak at the bond length."""
    frames = [
        Atoms(
            "Si2", positions=[(0, 0, 0), (2.35 + jitter, 0, 0)], cell=[10] * 3, pbc=True
        )
        for jitter in np_rng.normal(scale=0.01, size=20)
    ]
    radii, g_r = md_metrics.calc_rdf(frames, n_bins=100)
    assert radii[np.argmax(g_r)] == pytest.approx(2.35, abs=0.05)


def test_calc_rdf_r_max_respects_shrunken_cells() -> None:
    """Default r_max must come from the smallest cell; larger explicit ones raise."""
    frames = [
        Atoms("H2", positions=[(0, 0, 0), (1, 0, 0)], cell=[10] * 3, pbc=True),
        Atoms("H2", positions=[(0, 0, 0), (1, 0, 0)], cell=[4] * 3, pbc=True),
    ]

    radii, _g_r = md_metrics.calc_rdf(frames, n_bins=10)
    assert radii[-1] < 2  # default grid limited by the 4 A cell, not the 10 A one

    with pytest.raises(ValueError, match="minimum-image validity limit"):
        md_metrics.calc_rdf(frames, r_max=3)


def test_min_image_radius_uses_shortest_lattice_vector() -> None:
    """MIC radius is half the shortest lattice vector (Minkowski), not half the
    shortest basis vector, so skewed/non-reduced cells aren't over-cut and g(r) stays
    valid. A cell whose basis vectors are all ~5 A but whose b-a lattice vector is
    0.5 A must give radius 0.25 A, not 2.5 A.
    """
    pbc = np.array([True, True, True])
    cubic = (np.eye(3) * 5.0)[None]  # (1, 3, 3)
    assert md_metrics.min_image_radius(cubic, pbc) == pytest.approx(2.5)
    skewed = np.array([[[5.0, 0, 0], [5.0, 0.5, 0], [0, 0, 5.0]]])
    assert md_metrics.min_image_radius(skewed, pbc) == pytest.approx(0.25)

    # calc_rdf's default r_max + validation follow the corrected (smaller) limit
    frames = [
        Atoms("H2", positions=[[0, 0, 0], [0.05 * idx, 0, 0]], cell=skewed[0], pbc=True)
        for idx in range(4)
    ]
    radii, _ = md_metrics.calc_rdf(frames, n_bins=10)
    assert radii[-1] < 0.25  # grid capped at the shortest-lattice-vector limit
    with pytest.raises(ValueError, match="minimum-image validity limit"):
        md_metrics.calc_rdf(frames, n_bins=10, r_max=2.0)  # old basis-vector limit


@pytest.mark.parametrize(
    ("formulas", "n_bins", "err_msg"),
    [
        ([], 50, "zero frames"),
        (["H"], 50, "n_atoms=1 < 2"),
        (["H2", "H3"], 50, "Inconsistent atom counts"),
        # < 2 bins must raise instead of breaking downstream grid math
        (["H2"], 1, "Need >= 2 RDF bins"),
    ],
)
def test_calc_rdf_invalid_input(formulas: list[str], n_bins: int, err_msg: str) -> None:
    """RDF should reject bad bin counts and empty/single-atom/varying trajectories."""
    frames = [Atoms(formula, cell=[5] * 3, pbc=True) for formula in formulas]
    with pytest.raises(ValueError, match=err_msg):
        md_metrics.calc_rdf(frames, n_bins=n_bins)


@pytest.mark.parametrize(
    ("g_r_ref", "g_r_pred", "expected"),
    [
        ([0, 2, 1], [0, 2, 1], 0),  # perfect match
        ([0, 2, 1], [1, 1, 1], 100),  # prediction equals ideal gas
        ([1, 1, 1], [0, 2, 1], 100),  # reference equals ideal gas (zero denominator)
        # trapezoid(|3-5|)=2 vs trapezoid(|[0,3,0]-1|)=3 -> 2/3
        ([0, 3, 0], [0, 5, 0], 200 / 3),
        ([0, 1.5, 1], [50, 50, 50], 100),  # capped at 100%
    ],
)
def test_calc_rdf_error(
    g_r_ref: list[float], g_r_pred: list[float], expected: float
) -> None:
    """RDF error should match hand-computed values on small grids."""
    radii = np.array([1.0, 2.0, 3.0])
    err = md_metrics.calc_rdf_error(radii, np.array(g_r_ref), np.array(g_r_pred))
    assert err == pytest.approx(expected)

    with pytest.raises(ValueError, match="differ"):  # mismatched grid lengths
        md_metrics.calc_rdf_error(radii, np.ones(3), np.ones(4))


# === ADF ===


def right_angle_h3_frame(*, cell_len: float = 6) -> Atoms:
    """Three-atom frame with a 90-degree angle at atom 0."""
    positions = [(0, 0, 0), (1, 0, 0), (0, 1, 0)]
    return Atoms("H3", positions=positions, cell=[cell_len] * 3, pbc=True)


def tetrahedral_ch4_frame() -> Atoms:
    """Five-atom frame with a tetrahedral angle around atom 0."""
    neighbor_dirs = [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]
    positions = np.vstack([[0, 0, 0], neighbor_dirs]) / np.sqrt(3)
    return Atoms("CH4", positions=positions, cell=[8] * 3, pbc=True)


def wrapped_right_angle_h3_frame() -> Atoms:
    """Three-atom frame whose 90-degree angle crosses periodic boundaries."""
    positions = [(0.1, 0.1, 0.1), (3.1, 0.1, 0.1), (0.1, 3.1, 0.1)]
    return Atoms("H3", positions=positions, cell=[4] * 3, pbc=True)


@pytest.mark.parametrize(
    ("frame", "expected_angle", "abs_tol"),
    [
        (right_angle_h3_frame(), 90.5, 1e-12),
        (tetrahedral_ch4_frame(), np.degrees(np.arccos(-1 / 3)), 0.5),
        (wrapped_right_angle_h3_frame(), 90.5, 1e-12),
    ],
    ids=["right_angle", "tetrahedral", "wrapped_right_angle"],
)
def test_calc_adf_peak_position(
    frame: Atoms, expected_angle: float, abs_tol: float
) -> None:
    """ADF should recover known angles, including minimum-image neighbor vectors."""
    # bond_tolerance=2 so the artificial ~1 A H-H/C-H separations count as bonded
    angles, adf = md_metrics.calc_adf([frame], n_bins=180, bond_tolerance=2.0)
    assert angles[np.argmax(adf)] == pytest.approx(expected_angle, abs=abs_tol)
    np.testing.assert_allclose(
        np.sum(adf * (angles[1] - angles[0])), 1, rtol=1e-12, atol=1e-12
    )


def test_calc_adf_rigid_transform_invariant() -> None:
    """ADF should be invariant to rigid rotations and translations."""
    positions = [
        [0.0, 0.0, 0.0],
        [1.0, 0.2, -0.1],
        [-0.3, 1.1, 0.4],
        [0.5, -0.7, 1.2],
        [-1.0, -0.4, 0.6],
    ]
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    frame = Atoms("H5", positions=positions, cell=[50] * 3, pbc=False)
    transformed = Atoms(
        "H5", positions=positions @ rotation.T + [7, -3, 2], cell=[50] * 3, pbc=False
    )

    angles, adf = md_metrics.calc_adf([frame], n_bins=90, bond_tolerance=5.0)
    transformed_angles, transformed_adf = md_metrics.calc_adf(
        [transformed], n_bins=90, bond_tolerance=5.0
    )
    np.testing.assert_allclose(transformed_angles, angles, rtol=0, atol=0)
    np.testing.assert_allclose(transformed_adf, adf, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize(
    ("frames", "n_bins", "bond_tolerance", "err_msg"),
    [
        ([], 180, 1.25, "zero frames"),
        (
            [Atoms("H2", positions=[(0, 0, 0), (1, 0, 0)], cell=[6] * 3)],
            180,
            1.25,
            "n_atoms=2 < 3",
        ),
        ([right_angle_h3_frame()], 1, 1.25, "Need >= 2 ADF bins"),
        ([right_angle_h3_frame()], 180, 0.0, "must be positive"),
        # tiny tolerance leaves no covalent bonds -> no angle triplets
        ([right_angle_h3_frame()], 180, 0.1, "no bonded neighbor angle pairs"),
    ],
)
def test_calc_adf_invalid_input(
    frames: list[Atoms], n_bins: int, bond_tolerance: float, err_msg: str
) -> None:
    """ADF should reject invalid inputs and empty neighborhoods."""
    with pytest.raises(ValueError, match=err_msg):
        md_metrics.calc_adf(frames, n_bins=n_bins, bond_tolerance=bond_tolerance)


@pytest.mark.parametrize(
    ("angles", "adf_ref", "adf_pred", "expected"),
    [
        # identical -> W1 = 0
        ([0.0, 90.0, 180.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0),
        # half the mass shifted 0->90 deg: W1 = 45 deg, background W1 = 90 deg -> 50%
        ([0.0, 90.0, 180.0], [1.0, 0.0, 0.0], [0.5, 0.5, 0.0], 50),
        # prediction equals the sin background ([0,1,0]): W1 == background W1 -> 100%
        ([0.0, 90.0, 180.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 100),
        # mass shifted further than the background (0->180 deg): capped at 100%
        ([0.0, 90.0, 180.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 100),
    ],
    ids=["identical", "half_shift", "equals_background", "disjoint"],
)
def test_calc_adf_error(
    angles: list[float], adf_ref: list[float], adf_pred: list[float], expected: float
) -> None:
    """ADF error = W1(ref,pred) / W1(ref, sin-background), capped at 100%."""
    angles_arr, adf_ref_arr, adf_pred_arr = map(np.array, (angles, adf_ref, adf_pred))
    assert md_metrics.calc_adf_error(angles_arr, adf_ref_arr, adf_pred_arr) == (
        pytest.approx(expected)
    )


@pytest.mark.parametrize(
    ("angles", "adf_ref", "adf_pred", "err_msg"),
    [
        (np.array([0.0, 90.0, 180.0]), np.ones(3), np.ones(4), "differ"),
        (np.array([0.0, 90.0, 180.0]), np.zeros(3), np.ones(3), "non-positive area"),
        (
            np.array([0.0, 90.0, 180.0]),
            np.array([1.0, np.nan, 0.0]),
            np.ones(3),
            "non-finite intensities",
        ),
        (
            np.array([0.0, 90.0, 180.0]),
            np.array([1.0, -0.1, 0.0]),
            np.ones(3),
            "negative intensities",
        ),
        (np.array([0.0]), np.ones(1), np.ones(1), "angle grid"),
        (np.array([0.0, 2.0, 3.0]), np.ones(3), np.ones(3), "angle grid"),
    ],
    ids=["length_mismatch", "zero_area", "nan", "negative", "one_point", "irregular"],
)
def test_calc_adf_error_invalid_input(
    angles: np.ndarray, adf_ref: np.ndarray, adf_pred: np.ndarray, err_msg: str
) -> None:
    """ADF error should reject invalid grids and distributions."""
    with pytest.raises(ValueError, match=err_msg):
        md_metrics.calc_adf_error(angles, adf_ref, adf_pred)


# === velocities and vDOS ===


def test_calc_velocities_linear_motion() -> None:
    """Constant-velocity motion should give that exact velocity at every frame."""
    velocity = np.array([0.02, -0.01, 0.005])  # A/fs
    time_step_fs = 2.0
    frames = [
        Atoms("H", positions=[velocity * time_step_fs * step], cell=[100] * 3, pbc=True)
        for step in range(10)
    ]
    velocities = md_metrics.calc_velocities(frames, time_step_fs=time_step_fs)
    np.testing.assert_allclose(
        velocities, np.tile(velocity, (10, 1, 1)), rtol=1e-12, atol=1e-14
    )


def test_equipartition_temperature() -> None:
    """Equipartition T recovers the velocities' temperature and scales as 1/dt**2 (the
    property that catches a mislabeled dt_fs: e.g. dt 5x too large -> T 25x too low).
    """
    rng = np.random.default_rng(0)
    n_atoms, n_frames, dt = 32, 6, 1.5  # fs
    velocity = rng.normal(scale=0.01, size=(n_atoms, 3))  # A/fs, constant per atom
    base = rng.uniform(0, 10, size=(n_atoms, 3))
    frames = [
        Atoms("Cu32", positions=base + velocity * (step * dt), cell=[80] * 3, pbc=False)
        for step in range(n_frames)
    ]
    temp = md_metrics.equipartition_temperature(frames, time_step_fs=dt)
    # independent expectation from the same constant velocities (linear motion -> exact)
    factor = 1e10 * units.J / units.kg  # 1 amu*(A/fs)^2 in eV
    ke_ev = 0.5 * atomic_masses[29] * (velocity**2).sum() * factor
    expected = 2 * ke_ev / (3 * n_atoms * units.kB)
    assert temp == pytest.approx(expected, rel=1e-9)
    # halving the assumed dt quadruples the inferred temperature
    half_dt = md_metrics.equipartition_temperature(frames, time_step_fs=dt / 2)
    assert half_dt == pytest.approx(4 * temp, rel=1e-9)


def test_calc_velocities_unwraps_pbc_crossings() -> None:
    """Atoms crossing periodic boundaries must not produce velocity spikes."""
    time_step_fs, box_len, speed = 1.0, 5.0, 0.4  # A/fs, crosses boundary
    frames = [  # positions wrapped into the box by the % operator
        Atoms(
            "H",
            positions=[((speed * step) % box_len, 0, 0)],
            cell=[box_len] * 3,
            pbc=True,
        )
        for step in range(20)
    ]
    velocities = md_metrics.calc_velocities(frames, time_step_fs=time_step_fs)
    np.testing.assert_allclose(velocities[:, 0, 0], speed, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize(
    ("n_frames", "time_step_fs", "err_msg"),
    [(1, 1.0, "Need >= 2 frames"), (5, 0, "must be positive")],
)
def test_calc_velocities_invalid_input(
    n_frames: int, time_step_fs: float, err_msg: str
) -> None:
    """Velocity estimation should validate frame count and time step."""
    frames = [Atoms("H", positions=[(0, 0, 0)], cell=[5] * 3, pbc=True)] * n_frames
    with pytest.raises(ValueError, match=err_msg):
        md_metrics.calc_velocities(frames, time_step_fs=time_step_fs)


def make_cosine_velocities(
    *freq_amp_pairs: tuple[float, float],
    n_frames: int = 2000,
    time_step_fs: float = 2.0,
) -> np.ndarray:
    """Velocities of shape (n_frames, n_atoms, 3) where atom i oscillates along x
    with (frequency in THz, amplitude) = freq_amp_pairs[i].
    """
    times_ps = np.arange(n_frames) * time_step_fs * 1e-3
    velocities = np.zeros((n_frames, len(freq_amp_pairs), 3))
    for atom_idx, (freq_thz, amplitude) in enumerate(freq_amp_pairs):
        velocities[:, atom_idx, 0] = amplitude * np.cos(2 * np.pi * freq_thz * times_ps)
    return velocities


def test_calc_vdos_spectrum() -> None:
    """VDOS must peak at each atom's oscillation frequency, and per-DOF
    normalization (equipartition) must give a quiet atom's peak the same
    integrated weight as a 10x louder atom's peak.
    """
    velocities = make_cosine_velocities((5.0, 1.0), (12.5, 0.1))  # loud + quiet atom
    freqs, vdos = md_metrics.calc_vdos(velocities, time_step_fs=2.0)

    low_band = (freqs > 2.5) & (freqs < 7.5)
    high_band = (freqs > 10) & (freqs < 15)
    freq_resolution = freqs[1] - freqs[0]
    assert freqs[low_band][np.argmax(vdos[low_band])] == pytest.approx(
        5.0, abs=freq_resolution
    )
    assert freqs[high_band][np.argmax(vdos[high_band])] == pytest.approx(
        12.5, abs=freq_resolution
    )
    # a raw (unnormalized) power sum would give a weight ratio of ~100
    assert vdos[low_band].sum() / vdos[high_band].sum() == pytest.approx(1, rel=0.05)


@pytest.mark.parametrize(
    ("shape", "time_step_fs", "err_msg"),
    [
        ((10, 3), 1.0, "Expected shape"),
        # hanning(n < 4) windows are degenerate (all-zero for n=2), so 3 frames
        # must be rejected even though velocities only need 2
        ((3, 2, 3), 1.0, "Need >= 4 frames"),
        ((10, 2, 3), -1.0, "must be positive"),
    ],
)
def test_calc_vdos_invalid_input(
    shape: tuple[int, ...], time_step_fs: float, err_msg: str
) -> None:
    """VDOS should validate velocity array shape and time step."""
    with pytest.raises(ValueError, match=err_msg):
        md_metrics.calc_vdos(np.zeros(shape), time_step_fs=time_step_fs)


def test_calc_vdos_error_identical_and_disjoint() -> None:
    """Identical spectra give 0% error, non-overlapping spectra 100%."""
    freqs = np.linspace(0, 10, 101)
    peak_low = np.exp(-((freqs - 2) ** 2))
    peak_high = np.exp(-((freqs - 8) ** 2))

    assert md_metrics.calc_vdos_error(freqs, peak_low, freqs, peak_low) == 0
    err_disjoint = md_metrics.calc_vdos_error(freqs, peak_low, freqs, peak_high)
    assert err_disjoint == pytest.approx(100, abs=0.1)
    # partial overlap should land strictly in between
    err_partial = md_metrics.calc_vdos_error(
        freqs, peak_low, freqs, (peak_low + peak_high) / 2
    )
    assert 0 < err_partial < err_disjoint


def test_calc_vdos_error_interpolates_between_grids() -> None:
    """Same spectrum sampled on different grids should give near-zero error."""
    freqs_coarse = np.linspace(0, 10, 51)
    freqs_fine = np.linspace(0, 10, 201)
    spectrum = lambda freqs: np.exp(-((freqs - 5) ** 2))  # noqa: E731

    err = md_metrics.calc_vdos_error(
        freqs_coarse, spectrum(freqs_coarse), freqs_fine, spectrum(freqs_fine)
    )
    assert err == pytest.approx(0, abs=0.5)


@pytest.mark.parametrize(
    ("vdos_ref", "err_msg"),
    [(np.full(11, -1.0), "negative intensities"), (np.zeros(11), "non-positive area")],
)
def test_calc_vdos_error_invalid_spectra(vdos_ref: np.ndarray, err_msg: str) -> None:
    """VDOS error should reject negative or zero-area spectra."""
    freqs = np.linspace(0, 10, 11)
    with pytest.raises(ValueError, match=err_msg):
        md_metrics.calc_vdos_error(freqs, vdos_ref, freqs, np.ones(11))


# === pressure ===


def test_get_trajectory_pressures() -> None:
    """Per-frame pressures read from attached stress, averaging the full 3D trace (not a
    2D slab formula): the anisotropic frame gives 3 GPa, not 1.5.
    """
    pressures = [0.5, -1.0, 2.0]
    frames = [
        make_jiggled_frames(1, pressure_gpa=pressure)[0] for pressure in pressures
    ]
    aniso = Atoms("Si8", positions=np_rng.normal(scale=0.05, size=(8, 3)), cell=[6] * 3)
    aniso.pbc = True
    aniso.calc = SinglePointCalculator(aniso, stress=ANISOTROPIC_STRESS_VOIGT)
    np.testing.assert_allclose(
        md_metrics.get_trajectory_pressures([*frames, aniso]),
        [*pressures, 3.0],
        rtol=1e-10,
        atol=1e-12,
    )


@pytest.mark.parametrize(
    ("p_ref", "p_pred", "ref_dt", "pred_dt", "expected_mae"),
    [
        (np.array([0.0, 1.0, 2.0, 3.0]), np.array([0.0, 1.0, 2.0, 3.0]), 1, 1, 0),
        # constant +0.5 offset between equal-cadence trajectories
        (np.array([0.0, 1.0, 2.0, 3.0]), np.array([0.5, 1.5, 2.5, 3.5]), 1, 1, 0.5),
        # unequal lengths truncate the MAE to the shorter trajectory
        (np.array([0.0, 1.0, 2.0, 3.0]), np.array([0.5, 1.5]), 1, 1, 0.5),
        # differing cadences pair identical timestamps, not raw indices: pred samples
        # ref at every 2nd timestamp plus a pure +10 offset, so aligned MAE is exactly
        # 10 (index pairing would wrongly give 11.5)
        (np.arange(7.0), np.arange(7.0)[::2] + 10, 1, 2, 10),
        # awkward cadence ratio ref dt=0.3, pred dt=0.7 -> common grid every 2.1 fs;
        # common timestamps are 0 and 2.1 fs (ref indices 0,7; pred indices 0,3)
        (np.arange(8.0), np.array([100.0, -1.0, -2.0, 107.0]), 0.3, 0.7, 100),
    ],
    ids=["identical", "constant_offset", "truncated", "cadence_2x", "awkward_cadence"],
)
def test_calc_pressure_metrics_mae(
    p_ref: np.ndarray,
    p_pred: np.ndarray,
    ref_dt: float,
    pred_dt: float,
    expected_mae: float,
) -> None:
    """Pressure MAE pairs frames at identical timestamps, even across differing
    cadences and trajectory lengths.
    """
    metrics = md_metrics.calc_pressure_metrics(
        p_ref, p_pred, ref_time_step_fs=ref_dt, pred_time_step_fs=pred_dt
    )
    assert metrics["pressure_mae"] == pytest.approx(expected_mae)


def test_calc_pressure_metrics_wasserstein_and_validation() -> None:
    """W1 distance equals the mean offset and is frame-order-independent (unlike MAE);
    empty pressure arrays raise.
    """
    p_ref = np.array([0.0, 1.0, 2.0, 3.0])
    p_shifted = p_ref + 0.5

    assert md_metrics.calc_pressure_metrics(p_ref, p_ref)["pressure_wasserstein"] == 0
    shifted = md_metrics.calc_pressure_metrics(p_ref, p_shifted)
    assert shifted["pressure_wasserstein"] == pytest.approx(0.5)
    # W1 is invariant to frame order, MAE is not
    reversed_order = md_metrics.calc_pressure_metrics(p_ref, p_shifted[::-1])
    assert reversed_order["pressure_wasserstein"] == pytest.approx(0.5)
    assert reversed_order["pressure_mae"] > 0.5

    with pytest.raises(ValueError, match="empty"):
        md_metrics.calc_pressure_metrics(np.array([]), p_ref)


# === combined error ===


@pytest.mark.parametrize(
    ("rdf_error", "adf_error", "vdos_error", "pressure_error", "expected"),
    [
        (0, 0, 0, 0, 0),
        (10, 20, 30, 40, 25),  # simple mean of the four
        (15, 15, 15, 15, 15),
        (30, 0, 0, 0, 7.5),
    ],
)
def test_calc_combined_error(
    rdf_error: float,
    adf_error: float,
    vdos_error: float,
    pressure_error: float,
    expected: float,
) -> None:
    """Combined error: simple mean of RDF/ADF/vDOS/pressure errors."""
    assert md_metrics.calc_combined_error(
        rdf_error, adf_error, vdos_error, pressure_error
    ) == pytest.approx(expected)


def test_calc_combined_error_rejects_invalid_inputs() -> None:
    """Combined error should reject negative or non-finite components."""
    with pytest.raises(ValueError, match="must be non-negative"):
        md_metrics.calc_combined_error(-1, 5, 5, 5)
    # a missing/NaN pressure must fail loud, not silently become a two-metric mean
    with pytest.raises(ValueError, match="finite"):
        md_metrics.calc_combined_error(10, 20, 30, float("nan"))


@pytest.mark.parametrize(
    ("shift", "expected"),
    [(0, 0), (1000, 100)],  # identical -> 0% non-overlap; far apart -> ~100%
    ids=["identical", "disjoint"],
)
def test_calc_pressure_histogram_error(shift: float, expected: float) -> None:
    """E_P (Eq. 9): 0% for identical pressure samples, ~100% for disjoint supports."""
    pressures = np_rng.normal(0, 1, 5000)
    err = md_metrics.calc_pressure_histogram_error(pressures, pressures + shift)
    assert err == pytest.approx(expected, abs=1)


def test_calc_pressure_histogram_error_partial_and_symmetry() -> None:
    """Half-overlapping samples give 50% non-overlap; E_P is symmetric in its args and
    treats two identical single-valued distributions (hi == lo) as full overlap.
    """
    p_ref = np.array([0.25, 0.75, 1.25, 1.75])
    p_pred = np.array([1.25, 1.75, 2.25, 2.75])
    err = md_metrics.calc_pressure_histogram_error(p_ref, p_pred, n_bins=6)
    assert err == pytest.approx(50)
    assert md_metrics.calc_pressure_histogram_error(p_pred, p_ref, n_bins=6) == err
    const = np.full(8, 3.0)
    assert md_metrics.calc_pressure_histogram_error(const, const) == 0


def test_calc_pressure_histogram_error_rejects_bad_input() -> None:
    """Empty arrays, <2 bins and non-finite pressures raise, not return garbage."""
    pressures = np_rng.normal(0, 1, 100)
    with pytest.raises(ValueError, match="empty"):
        md_metrics.calc_pressure_histogram_error(np.array([]), pressures)
    with pytest.raises(ValueError, match="bins"):
        md_metrics.calc_pressure_histogram_error(pressures, pressures, n_bins=1)
    with pytest.raises(ValueError, match="finite"):
        md_metrics.calc_pressure_histogram_error(np.array([0.0, np.nan]), pressures)


# === energy/force RMSE ===


@pytest.mark.parametrize("force_offset", [0, 0.1, 0.25])
def test_calc_energy_force_rmse(force_offset: float) -> None:
    """Energy RMSE uses mean-subtracted (fluctuation) energies, so a constant predictor
    vs constant references cancels to 0 regardless of the gap; forces map directly to
    the absolute error (they're invariant to the energy zero).
    """
    frames = [h2_frame(0.0) for _ in range(3)]
    calc = ConstantCalculator(energy=0.8, force=force_offset)
    result = md_metrics.calc_energy_force_rmse(frames, calc)
    assert result["energy_rmse"] == pytest.approx(0, abs=1e-12)  # constant gap removed
    assert result["force_rmse"] == pytest.approx(abs(force_offset), abs=1e-12)

    with pytest.raises(ValueError, match="zero frames"):
        md_metrics.calc_energy_force_rmse([], calc)
    single_frame = md_metrics.calc_energy_force_rmse(frames[:1], calc)
    assert single_frame["energy_rmse"] == pytest.approx(0, abs=1e-12)  # no fluctuation


def test_calc_energy_force_rmse_removes_offset_keeps_fluctuations() -> None:
    """Energy RMSE subtracts each trajectory's mean, so a large constant reference
    offset drops out and only the per-frame fluctuation mismatch remains. This is what
    makes it invariant to CFPMD-26's mixed all-electron/PAW absolute energy zeros.
    """
    n_atoms = 2
    fluctuations = np.array([0.0, 0.3, -0.3, 0.6, -0.6])  # eV total, varying, zero mean
    # constant predictor returns 0; references carry a +1000 eV offset plus fluctuations
    frames = [h2_frame(-(1000 + float(fluct))) for fluct in fluctuations]

    result = md_metrics.calc_energy_force_rmse(frames, ConstantCalculator(energy=0))
    # the 1000 eV offset cancels; the flat predictor's deviation is 0 and the reference
    # deviation is -fluct (mean(fluctuations)=0), so only the fluctuations survive
    expected = np.sqrt(np.mean((fluctuations / n_atoms) ** 2))
    assert result["energy_rmse"] == pytest.approx(expected)


def test_calc_energy_force_rmse_uses_frame_geometry() -> None:
    """The calculator must receive each frame's own geometry (not a shared/blank one),
    so a geometry-dependent prediction yields per-frame-varying errors.
    """

    class XCoordCalculator(Calculator):
        """Energy equals atom 0's x-coordinate, so it depends on the frame geometry."""

        implemented_properties = ("energy", "forces")

        def calculate(
            self,
            atoms: Atoms | None = None,
            properties: list[str] | None = None,
            system_changes: list[str] = all_changes,
        ) -> None:
            super().calculate(atoms, properties, system_changes)
            assert atoms is not None
            self.results["energy"] = float(atoms.positions[0, 0])
            self.results["forces"] = np.zeros((len(atoms), 3))

    x_coords = np.array([0.0, 1.0, 2.0, 3.0])
    frames = [
        h2_frame(0.0, positions=np.array([[x, 0, 0], [0, 0, 0]])) for x in x_coords
    ]
    result = md_metrics.calc_energy_force_rmse(frames, XCoordCalculator())
    # pred energy = x (ref is flat 0), so after mean-subtraction the per-frame deviation
    # is (x - mean(x)) / 2 atoms
    expected = np.sqrt(np.mean(((x_coords - x_coords.mean()) / 2) ** 2))
    assert result["energy_rmse"] == pytest.approx(expected)


def test_predict_on_reference_persists_for_cpu_recompute(tmp_path: Path) -> None:
    """Persisted per-frame predictions recompute the same energy/force RMSE on CPU as
    the calc_energy_force_rmse wrapper, with no model re-run.
    """
    frames = [
        h2_frame(float(np_rng.random()), forces=np_rng.random((2, 3))) for _ in range(6)
    ]
    calc = ConstantCalculator(energy=0.5, force=0.1)
    preds = md_metrics.predict_on_reference(frames, calc)
    assert set(preds) == {"e_pred", "force_se"}
    assert preds["e_pred"].shape == preds["force_se"].shape == (6,)

    direct = md_metrics.calc_energy_force_rmse(frames, calc)
    assert md_metrics.energy_force_rmse_from_preds(frames, preds) == direct

    # persisted sidecar recomputes identical metrics without the calculator
    path = tmp_path / "refeval.npz"
    np.savez(path, e_pred=preds["e_pred"], force_se=preds["force_se"])
    with np.load(path) as data:
        reloaded = {key: data[key] for key in data.files}
    assert md_metrics.energy_force_rmse_from_preds(frames, reloaded) == direct

    # a stale sidecar with the wrong frame count must error, not silently miscompute
    stale = {key: arr[:-1] for key, arr in preds.items()}
    with pytest.raises(ValueError, match="stale or mismatched"):
        md_metrics.energy_force_rmse_from_preds(frames, stale)


def test_calc_rdf_matches_direct_ase_histogram() -> None:
    """The array-based calc_rdf reproduces a direct Atoms.get_all_distances(mic=True)
    histogram for a triclinic cell (independent cross-check of the new kernel).
    """
    cell = np.array([[6.0, 0, 0], [0.8, 6.0, 0], [0.3, -0.4, 6.0]])  # triclinic
    frames = [
        Atoms("Cu8", positions=np_rng.random((8, 3)) * 5, cell=cell, pbc=True)
        for _ in range(3)
    ]
    n_bins = 40
    radii, g_r = md_metrics.calc_rdf(frames, n_bins=n_bins)

    r_max = float(radii[-1] + (radii[1] - radii[0]) / 2)
    bin_edges = np.linspace(0, r_max, n_bins + 1)
    counts = np.zeros(n_bins)
    inv_volume_sum = 0.0
    upper_tri = np.triu_indices(8, k=1)
    for atoms in frames:
        counts += np.histogram(
            atoms.get_all_distances(mic=True)[upper_tri], bins=bin_edges
        )[0]
        inv_volume_sum += 1 / atoms.get_volume()
    shell_volumes = 4 / 3 * np.pi * np.diff(bin_edges**3)
    g_r_direct = counts / (8**2 / 2 * inv_volume_sum * shell_volumes)
    np.testing.assert_allclose(g_r, g_r_direct, rtol=0, atol=1e-12)


# === frame matching and per-system evaluation ===


@pytest.mark.parametrize(
    ("n_ref", "n_pred", "ref_dt", "pred_dt", "expected"),
    [
        # spans are (n_frames - 1) * dt; both outputs cover identical [0, t_match]
        (10, 10, 1, 1, (10, 10)),  # equal everything
        (10, 20, 2, 1, (10, 19)),  # ref spans 18 fs < pred's 19 fs
        (10, 5, 1, 1, (5, 5)),  # pred shorter
        (4, 100, 2.5, 1, (3, 6)),  # ref spans 7.5 fs, snaps to 5 fs common grid
        (3, 100, 2.5, 1, (3, 6)),  # ref spans exactly 5 fs = one common step
        (7, 100, 2.5, 1, (7, 16)),  # ref spans exactly 15 fs = three common steps
        (1, 10, 1, 1, (1, 1)),  # single frame matches only t=0
    ],
)
def test_matched_frame_counts(
    n_ref: int, n_pred: int, ref_dt: float, pred_dt: float, expected: tuple[int, int]
) -> None:
    """Matched frame counts should span the same total simulation time."""
    result = md_metrics.matched_frame_counts(
        n_ref_frames=n_ref,
        n_pred_frames=n_pred,
        ref_time_step_fs=ref_dt,
        pred_time_step_fs=pred_dt,
    )
    assert result == expected

    with pytest.raises(ValueError, match="must be positive"):  # zero time step
        md_metrics.matched_frame_counts(
            n_ref_frames=5, n_pred_frames=5, ref_time_step_fs=0, pred_time_step_fs=1
        )


def test_evaluate_md_system() -> None:
    """Per-system evaluation should produce all metrics on matched trajectories."""
    ref_traj = make_jiggled_frames(10, pressure_gpa=1.0)
    pred_traj = make_jiggled_frames(20, pressure_gpa=1.2)

    metrics = md_metrics.evaluate_md_system(
        ref_traj,
        pred_traj,
        ref_time_step_fs=2,
        pred_time_step_fs=1,
        n_rdf_bins=50,
    )

    assert set(metrics) == {
        "rdf_error",
        "adf_error",
        "vdos_error",
        "pressure_mae",
        "pressure_wasserstein",
        "pressure_error",
    }
    assert 0 <= metrics["rdf_error"] <= 100
    assert 0 <= metrics["adf_error"] <= 100
    assert 0 <= metrics["vdos_error"] <= 100
    assert metrics["pressure_mae"] == pytest.approx(0.2, abs=1e-10)
    assert metrics["pressure_wasserstein"] == pytest.approx(0.2, abs=1e-10)
    assert 0 <= metrics["pressure_error"] <= 100

    # frames without stress data should yield NaN pressure metrics, not errors
    no_stress = md_metrics.evaluate_md_system(
        make_jiggled_frames(8),
        make_jiggled_frames(8),
        ref_time_step_fs=1,
        pred_time_step_fs=1,
    )
    assert np.isnan(no_stress["pressure_mae"])
    assert np.isnan(no_stress["pressure_wasserstein"])
    assert np.isnan(no_stress["pressure_error"])


def test_evaluate_md_system_reraises_unexpected_pressure_errors() -> None:
    """Unexpected stress backend failures should not become NaN pressure metrics."""
    ref_traj = make_jiggled_frames(8)
    pred_traj = make_jiggled_frames(8)
    for atoms in (*ref_traj, *pred_traj):
        atoms.calc = BrokenStressCalculator()

    with pytest.raises(RuntimeError, match="stress backend failure"):
        md_metrics.evaluate_md_system(
            ref_traj,
            pred_traj,
            ref_time_step_fs=1,
            pred_time_step_fs=1,
        )


def test_evaluate_md_system_with_calculator_adds_rmse() -> None:
    """Passing a calculator should add energy and force RMSE metrics."""
    ref_traj = [
        attach_ref_labels(atoms, energy=-1.0, forces=np.zeros((len(atoms), 3)))
        for atoms in make_jiggled_frames(6)
    ]

    metrics = md_metrics.evaluate_md_system(
        ref_traj,
        make_jiggled_frames(6),
        ref_time_step_fs=1,
        pred_time_step_fs=1,
        calculator=ConstantCalculator(energy=-0.2),  # constant gap vs flat ref -1.0
    )
    assert metrics["energy_rmse"] == pytest.approx(0)  # constant offset cancels
    assert metrics["force_rmse"] == pytest.approx(0)

    # trajectories shorter than two matched frames should raise
    with pytest.raises(ValueError, match="too short after time matching"):
        md_metrics.evaluate_md_system(
            make_jiggled_frames(1),
            make_jiggled_frames(10),
            ref_time_step_fs=1,
            pred_time_step_fs=1,
        )


# === aggregation and YAML output ===


def test_calc_md_metrics() -> None:
    """Aggregation should average per-system rows, skipping NaNs."""
    df_md = pd.DataFrame(
        {
            "energy_rmse": [0.1, 0.3],
            "force_rmse": [0.2, 0.4],
            "rdf_error": [10.0, 20.0],
            "adf_error": [50.0, 10.0],
            "vdos_error": [30.0, 10.0],
            "pressure_mae": [1.0, np.nan],  # one system without stress data
            "pressure_wasserstein": [0.5, np.nan],
            "pressure_error": [40.0, np.nan],
        }
    )
    metrics = md_metrics.calc_md_metrics(df_md)

    # energy/force RMSE means (0.2, 0.3 eV) are reported in meV (x1000)
    assert metrics["energy_rmse"] == pytest.approx(200)
    assert metrics["force_rmse"] == pytest.approx(300)
    assert metrics["rdf_error"] == pytest.approx(15)
    assert metrics["adf_error"] == pytest.approx(30)
    assert metrics["vdos_error"] == pytest.approx(20)
    assert metrics["pressure_mae"] == pytest.approx(1)
    assert metrics["pressure_wasserstein"] == pytest.approx(0.5)
    assert metrics["pressure_error"] == pytest.approx(40)
    # combined = simple mean of mean rdf/adf/vdos/pressure errors
    assert metrics["combined_error"] == pytest.approx(
        md_metrics.calc_combined_error(15, 30, 20, 40)
    )
    assert metrics["n_systems"] == 2

    # partial columns work, unrecognized columns raise
    partial = md_metrics.calc_md_metrics(pd.DataFrame({"rdf_error": [10.0, 20.0]}))
    assert partial == {"rdf_error": 15, "n_systems": 2}
    with pytest.raises(ValueError, match="No recognized MD metric columns"):
        md_metrics.calc_md_metrics(pd.DataFrame({"unrelated": [1]}))


def test_calc_md_metrics_skips_combined_error_without_finite_pressure() -> None:
    """All-NaN pressure errors should not make aggregation fail or emit CMDS."""
    df_md = pd.DataFrame(
        {
            "rdf_error": [10.0, 20.0],
            "adf_error": [30.0, 40.0],
            "vdos_error": [50.0, 60.0],
            "pressure_error": [np.nan, np.nan],
        }
    )

    metrics = md_metrics.calc_md_metrics(df_md)

    assert metrics["rdf_error"] == pytest.approx(15)
    assert metrics["adf_error"] == pytest.approx(35)
    assert metrics["vdos_error"] == pytest.approx(55)
    assert np.isnan(metrics["pressure_error"])
    assert "combined_error" not in metrics


def test_combine_per_system_metrics() -> None:
    """Per-system frames concat into one frame indexed by system; reruns override."""
    frames = [
        pd.DataFrame({"system": ["sysA"], "rdf_error": [10.0]}),
        pd.DataFrame({"system": ["sysB"], "rdf_error": [20.0]}),
        # a rerun of sysA with a better value, later in the list -> must override
        pd.DataFrame({"system": ["sysA"], "rdf_error": [5.0]}),
    ]
    df_md = md_metrics.combine_per_system_metrics(frames)
    assert list(df_md.index) == ["sysB", "sysA"]  # sysA dedup'd to its last (rerun) row
    assert df_md.loc["sysA", "rdf_error"] == 5.0
    assert df_md.loc["sysB", "rdf_error"] == 20.0

    with pytest.raises(ValueError, match="No per-system metric frames given"):
        md_metrics.combine_per_system_metrics([])
    with pytest.raises(ValueError, match="lack a 'system' column"):
        md_metrics.combine_per_system_metrics([pd.DataFrame({"rdf_error": [10.0]})])


@pytest.mark.parametrize(
    "init_yaml",  # placeholder strings at leaf and intermediate level get replaced
    ["metrics:\n  md: not available\n", "metrics: not available\n"],
)
def test_write_metrics_to_yaml(
    tmp_path: Path, init_yaml: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Metrics should be written under metrics.md with units and rounding."""
    yaml_path = tmp_path / "model.yml"
    yaml_path.write_text(f"model_name: test\n{init_yaml}", encoding="utf-8")
    model = Model.mace_mp_0
    monkeypatch.setattr(Model, "yaml_path", yaml_path)

    metrics = {
        "energy_rmse": 0.123456,
        "force_rmse": 0.234567,
        "rdf_error": 12.34567,
        "combined_error": 23.456789,
        "n_systems": 17,
    }
    path = "models/test/md-metrics.csv"
    url = "https://example.com/md-metrics.csv"
    md_metrics.write_metrics_to_yaml(
        model, metrics, pred_file_path=path, pred_file_url=url
    )

    # read back as UTF-8 so the Å in the force_rmse unit decodes correctly on Windows
    text = yaml_path.read_text(encoding="utf-8")
    assert "pred_file: models/test/md-metrics.csv" in text
    assert "pred_file_url: https://example.com/md-metrics.csv" in text
    assert f"energy_rmse: 0.1235 # {md_metrics.METRIC_UNITS['energy_rmse']}" in text
    assert f"force_rmse: 0.2346 # {md_metrics.METRIC_UNITS['force_rmse']}" in text
    assert "rdf_error: 12.3457 # %" in text
    assert "combined_error: 23.4568 # %" in text
    assert "n_systems: 17 # count" in text
    # the 'not available' placeholder must be replaced, not kept
    assert "not available" not in text
