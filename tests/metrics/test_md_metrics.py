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
    """Predicts a fixed energy and uniform force regardless of geometry."""

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
    """An H2 reference frame carrying private energy/force labels."""
    if positions is None:
        positions = np_rng.random((2, 3)) * 2
    if forces is None:
        forces = np.zeros((2, 3))
    atoms = Atoms("H2", positions=positions, cell=[5] * 3)
    atoms.calc = SinglePointCalculator(atoms, energy=energy, forces=forces)
    return atoms


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
    with pytest.raises(ValueError, match="negative"):
        md_metrics.calc_rdf_error(radii, np.array([0.0, -1.0, 1.0]), np.ones(3))
    with pytest.raises(ValueError, match="non-finite"):
        md_metrics.calc_rdf_error(radii, np.array([0.0, np.inf, 1.0]), np.ones(3))


# === ADF ===


def right_angle_h3_frame(*, cell_len: float = 6, pbc: bool = True) -> Atoms:
    """Three-atom frame with a 90-degree angle at atom 0."""
    positions = [(0, 0, 0), (1, 0, 0), (0, 1, 0)]
    return Atoms("H3", positions=positions, cell=[cell_len] * 3, pbc=pbc)


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
        (right_angle_h3_frame(cell_len=1, pbc=False), 90.5, 1e-12),
    ],
    ids=["right", "tetrahedral", "wrapped", "non_pbc_small_cell"],
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
    assert md_metrics.calc_adf_error(
        np.array(angles), np.array(adf_ref), np.array(adf_pred)
    ) == pytest.approx(expected)


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


def test_calc_vdos_error_w1_shift_monotonic() -> None:
    """W1-based vDOS error: 0 for identical spectra, monotonic in peak shift (for equal
    Gaussians error% == 100*min(1, shift/sigma_ref)), capped at 100% beyond sigma_ref.
    """
    freqs = np.linspace(0, 40, 801)
    ref = np.exp(-((freqs - 10) ** 2) / 2)  # sigma_ref = 1 THz, centered at 10 THz
    shifted = lambda delta: np.exp(-((freqs - (10 + delta)) ** 2) / 2)  # noqa: E731

    assert md_metrics.calc_vdos_error(
        freqs, ref, freqs, ref, band_frac=1.0
    ) == pytest.approx(0, abs=1e-9)
    errs = [
        md_metrics.calc_vdos_error(freqs, ref, freqs, shifted(delta), band_frac=1.0)
        for delta in (0.25, 0.5, 0.75)
    ]
    assert errs == pytest.approx([25, 50, 75], abs=1)  # error% == 100*shift/sigma_ref
    # a shift beyond the reference's own spectral width saturates at the 100% cap
    assert md_metrics.calc_vdos_error(
        freqs, ref, freqs, shifted(1.5), band_frac=1.0
    ) == pytest.approx(100)


@pytest.mark.parametrize(
    ("vdos_ref", "err_msg"),
    [
        (np.full(11, -1.0), "negative intensities"),
        (np.zeros(11), "non-positive area"),
        (np.full(11, np.inf), "non-finite intensities"),
    ],
)
def test_calc_vdos_error_invalid_spectra(vdos_ref: np.ndarray, err_msg: str) -> None:
    """Reject negative, zero-area or non-finite vDOS spectra."""
    freqs = np.linspace(0, 10, 11)
    with pytest.raises(ValueError, match=err_msg):
        md_metrics.calc_vdos_error(freqs, vdos_ref, freqs, np.ones(11))


@pytest.mark.parametrize("band_frac", [0.0, -0.1, 1.5])
def test_calc_vdos_error_invalid_band_frac(band_frac: float) -> None:
    """Reject band_frac outside (0, 1]."""
    freqs, ones = np.linspace(0, 10, 11), np.ones(11)
    with pytest.raises(ValueError, match="band_frac"):
        md_metrics.calc_vdos_error(freqs, ones, freqs, ones, band_frac=band_frac)


def test_calc_vdos_error_edge_cases() -> None:
    """A prediction with no power in the reference band caps at 100%; a degenerate
    single-frequency reference returns 0 only when the prediction is the same line.
    """
    freqs = np.linspace(0, 40, 401)
    ref = np.exp(-((freqs - 5) ** 2) / 2)  # low-frequency reference band
    pred_far = np.exp(-((freqs - 35) ** 2) / 2)  # all power above the clipped band
    assert md_metrics.calc_vdos_error(freqs, ref, freqs, pred_far) == pytest.approx(100)

    delta = np.zeros(401)
    delta[100] = 1.0  # single non-zero frequency bin -> sigma_ref == 0
    assert md_metrics.calc_vdos_error(freqs, delta, freqs, delta, band_frac=1.0) == 0
    other = np.zeros(401)
    other[200] = 1.0
    assert md_metrics.calc_vdos_error(
        freqs, delta, freqs, other, band_frac=1.0
    ) == pytest.approx(100)


def test_calc_vdos_error_band_clip_suppresses_tail() -> None:
    """Clipping (band_frac < 1) drops the reference's noisy tail, so a spurious far
    prediction peak there is penalized less than at band_frac=1.0 (full grid).
    """
    freqs = np.linspace(0, 40, 401)
    ref = np.exp(-((freqs - 8) ** 2) / 2)
    pred = ref + 0.3 * np.exp(-((freqs - 32) ** 2) / 2)  # spurious high-frequency peak
    clipped = md_metrics.calc_vdos_error(freqs, ref, freqs, pred, band_frac=0.999)
    unclipped = md_metrics.calc_vdos_error(freqs, ref, freqs, pred, band_frac=1.0)
    assert clipped < unclipped


def test_calc_vdos_error_grid_sampling_invariant() -> None:
    """The error reflects spectral mass, not merged-grid sampling: the same physical
    spectrum on different (incommensurate) rfftfreq grids scores ~0, and a shifted
    spectrum gives the same error whichever grid each side uses (raw-density weighting
    would drift with sampling density).
    """
    shape = lambda freq, mu: np.exp(-((freq - mu) ** 2) / 8)  # noqa: E731
    coarse = np.fft.rfftfreq(2001, 0.0025)  # ~5 ps trajectory, 2.5 fs step
    fine = np.fft.rfftfreq(5001, 0.001)  # ~5 ps trajectory, 1 fs step
    # identical physical spectrum sampled on the two grids -> near-zero error
    assert md_metrics.calc_vdos_error(
        coarse, shape(coarse, 10), fine, shape(fine, 10)
    ) == pytest.approx(0, abs=0.1)
    # a 1 THz peak shift gives the same error regardless of each side's grid
    converged = md_metrics.calc_vdos_error(fine, shape(fine, 10), fine, shape(fine, 11))
    mixed_a = md_metrics.calc_vdos_error(
        coarse, shape(coarse, 10), fine, shape(fine, 11)
    )
    mixed_b = md_metrics.calc_vdos_error(
        fine, shape(fine, 10), coarse, shape(coarse, 11)
    )
    assert mixed_a == pytest.approx(converged, abs=0.1)
    assert mixed_b == pytest.approx(converged, abs=0.1)


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
    ("p_ref", "p_pred", "expected_mae"),
    [
        (np.array([0.0, 1.0, 2.0, 3.0]), np.array([0.0, 1.0, 2.0, 3.0]), 0.0),
        # constant +0.5 offset -> mean bias 0.5
        (np.array([0.0, 1.0, 2.0, 3.0]), np.array([0.5, 1.5, 2.5, 3.5]), 0.5),
        # unequal lengths: means are 1.5 vs 1.0 -> bias 0.5
        (np.array([0.0, 1.0, 2.0, 3.0]), np.array([0.5, 1.5]), 0.5),
        # frame order is irrelevant: reversal preserves the mean -> zero bias
        (np.arange(7.0), np.arange(7.0)[::-1], 0.0),
    ],
    ids=["identical", "offset", "unequal_len", "reversed"],
)
def test_calc_pressure_metrics_mae(
    p_ref: np.ndarray, p_pred: np.ndarray, expected_mae: float
) -> None:
    """Pressure MAE is the mean-pressure bias |mean(ref) - mean(pred)|: it compares the
    separately-averaged trajectory pressures, independent of frame order and length.
    """
    metrics = md_metrics.calc_pressure_metrics(p_ref, p_pred)
    assert metrics["pressure_mae"] == pytest.approx(expected_mae)


def test_calc_pressure_metrics_wasserstein_and_validation() -> None:
    """W1 distance equals the mean offset; both W1 and the mean-bias MAE are frame-order
    independent (they compare averaged/full distributions); empty pressure arrays raise.
    """
    p_ref = np.array([0.0, 1.0, 2.0, 3.0])
    p_shifted = p_ref + 0.5

    assert md_metrics.calc_pressure_metrics(p_ref, p_ref)["pressure_wasserstein"] == 0
    shifted = md_metrics.calc_pressure_metrics(p_ref, p_shifted)
    assert shifted["pressure_wasserstein"] == pytest.approx(0.5)
    assert shifted["pressure_mae"] == pytest.approx(0.5)
    # both metrics are invariant to frame order (mean/distribution based)
    reversed_order = md_metrics.calc_pressure_metrics(p_ref, p_shifted[::-1])
    assert reversed_order["pressure_wasserstein"] == pytest.approx(0.5)
    assert reversed_order["pressure_mae"] == pytest.approx(0.5)

    with pytest.raises(ValueError, match="empty"):
        md_metrics.calc_pressure_metrics(np.array([]), p_ref)


# === private energy/force diagnostics ===


def test_calc_energy_force_rmse_uses_private_labels() -> None:
    """Energy fluctuations and force RMSE come from private reference labels."""
    frames = [
        h2_frame(energy=1.0, forces=np.zeros((2, 3))),
        h2_frame(energy=3.0, forces=np.ones((2, 3))),
    ]

    metrics = md_metrics.calc_energy_force_rmse(
        frames, ConstantCalculator(energy=2.0, force=0.0)
    )

    # predicted energy fluctuations are zero; reference fluctuations are ±0.5 eV/atom
    assert metrics["energy_rmse"] == pytest.approx(0.5)
    assert metrics["force_rmse"] == pytest.approx(np.sqrt(0.5))


def test_energy_force_rmse_rejects_missing_private_labels() -> None:
    """Private diagnostics must fail loudly when the reference lacks labels."""
    atoms = Atoms("H2", positions=np_rng.random((2, 3)) * 2, cell=[5] * 3)

    with pytest.raises(ValueError, match="energy/forces"):
        md_metrics.calc_energy_force_rmse([atoms], ConstantCalculator())


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


def test_evaluate_md_system_rejects_short_trajectories() -> None:
    """Trajectories shorter than two matched frames should raise."""
    with pytest.raises(ValueError, match="too short after time matching"):
        md_metrics.evaluate_md_system(
            make_jiggled_frames(1),
            make_jiggled_frames(10),
            ref_time_step_fs=1,
            pred_time_step_fs=1,
        )


# === aggregation and YAML output ===


def test_calc_md_metrics() -> None:
    """Aggregation should average per-system rows, skipping NaNs, sum run_time_sec,
    max the memory peaks and pass through a unique hardware label.
    """
    df_md = pd.DataFrame(
        {
            "energy_rmse": [0.001, 0.003],
            "force_rmse": [0.1, 0.3],
            "rdf_error": [10.0, 20.0],
            "adf_error": [50.0, 10.0],
            "vdos_error": [30.0, 10.0],
            "pressure_mae": [1.0, np.nan],  # one system without stress data
            "pressure_wasserstein": [0.5, np.nan],
            "pressure_error": [40.0, np.nan],
            "run_time_sec": [100.5, 200.25],
            "max_rss_gb": [4.2, 6.1],
            "max_gpu_mem_gb": [11.5, 9.0],
            "hardware": ["NVIDIA H200", "NVIDIA H200"],
        }
    )
    metrics = md_metrics.calc_md_metrics(df_md)

    assert metrics.pop("hardware") == "NVIDIA H200"
    # exact-key equality also proves combined_score is absent (CMDS is site-computed)
    assert metrics == pytest.approx(
        {
            "run_time_sec": 300.75,  # summed, not averaged
            "max_rss_gb": 6.1,  # max over systems, not averaged
            "max_gpu_mem_gb": 11.5,
            "energy_rmse": 2,
            "force_rmse": 200,
            "rdf_error": 15,
            "adf_error": 30,
            "vdos_error": 20,
            "pressure_mae": 1,
            "pressure_wasserstein": 0.5,
            "pressure_error": 40,
            "n_systems": 2,
        }
    )

    # partial columns work, unrecognized columns raise
    partial = md_metrics.calc_md_metrics(pd.DataFrame({"rdf_error": [10.0, 20.0]}))
    assert partial == {"rdf_error": 15, "n_systems": 2}
    with pytest.raises(ValueError, match="No recognized MD metric columns"):
        md_metrics.calc_md_metrics(pd.DataFrame({"unrelated": [1]}))


@pytest.mark.parametrize(
    ("cols", "dropped"),
    [
        (
            {"energy_rmse": [0.001, np.nan], "force_rmse": [0.1, np.nan]},
            ("energy_rmse", "force_rmse"),
        ),
        (
            {"run_time_sec": [100.0, np.nan], "hardware": ["NVIDIA H200"] * 2},
            ("run_time_sec",),
        ),
        (
            # mixed GPUs also drop the GPU-dependent costs: a wall-time sum or VRAM
            # max over an H200 and an A100 misrepresents the model's cost with no
            # hardware label left to qualify it. Host RSS is GPU-independent -> kept
            {
                "run_time_sec": [100.0, 200.0],
                "max_rss_gb": [4.2, 6.1],
                "max_gpu_mem_gb": [11.5, 9.0],
                "hardware": ["NVIDIA H200", "NVIDIA A100"],
            },
            ("hardware", "run_time_sec", "max_gpu_mem_gb"),
        ),
        (
            {"run_time_sec": [100.0, 200.0], "hardware": ["NVIDIA H200", None]},
            ("hardware", "run_time_sec"),
        ),
        (
            # partial memory coverage -> dropped (a subset max could hide the true
            # peak); the complete sibling memory column is kept
            {"max_rss_gb": [4.2, np.nan], "max_gpu_mem_gb": [11.5, 9.0]},
            ("max_rss_gb",),
        ),
    ],
    ids=[
        "partial_energy_force",
        "partial_run_time",
        "mixed_gpus",
        "missing_gpu",
        "partial_memory",
    ],
)
def test_calc_md_metrics_all_or_nothing_columns(
    cols: dict[str, list[float | str | None]], dropped: tuple[str, ...]
) -> None:
    """E/F diagnostics and run provenance publish only when every system row agrees:
    partial coverage or a mixed/missing GPU label is dropped (a subset mean/sum would
    misrepresent the model while n_systems suggests full coverage), but complete
    sibling columns are kept.
    """
    df_md = pd.DataFrame({"rdf_error": [10.0, 20.0]} | cols)
    metrics = md_metrics.calc_md_metrics(df_md)
    assert set(dropped).isdisjoint(metrics)
    assert set(cols) - set(dropped) <= set(metrics)
    assert metrics["rdf_error"] == pytest.approx(15)
    assert metrics["n_systems"] == 2


def test_calc_md_metrics_handles_all_nan_pressure() -> None:
    """All-NaN pressure errors should not make aggregation fail."""
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

    metrics: dict[str, float | str] = {
        "hardware": "NVIDIA H200",
        "run_time_sec": 4521.484,
        "energy_rmse": 1.23456,
        "force_rmse": 98.7654,
        "rdf_error": 12.34567,
        "n_systems": 17,
    }
    path = "models/test/md-metrics.csv"
    url = "https://example.com/md-metrics.csv"
    md_metrics.write_metrics_to_yaml(
        model, metrics, pred_file_path=path, pred_file_url=url
    )

    text = yaml_path.read_text(encoding="utf-8")
    assert "pred_file: models/test/md-metrics.csv" in text
    assert "pred_file_url: https://example.com/md-metrics.csv" in text
    assert "hardware: NVIDIA H200" in text
    assert "run_time_sec: 4521.484 # s" in text
    assert "energy_rmse: 1.2346 # meV/atom" in text
    assert "force_rmse: 98.7654 # meV/Å" in text
    assert "rdf_error: 12.3457 # %" in text
    assert "n_systems: 17 # count" in text
    # CMDS is site-computed, never persisted
    assert "combined_score" not in text
    # the 'not available' placeholder must be replaced, not kept
    assert "not available" not in text

    # a later recompute from per-system CSVs lacking timing columns (e.g. legacy runs)
    # must preserve the recorded run provenance, not silently drop it
    md_metrics.write_metrics_to_yaml(model, {"rdf_error": 11.0, "n_systems": 17})
    text = yaml_path.read_text(encoding="utf-8")
    assert "hardware: NVIDIA H200" in text
    assert "run_time_sec: 4521.484" in text
    assert "rdf_error: 11.0 # %" in text
