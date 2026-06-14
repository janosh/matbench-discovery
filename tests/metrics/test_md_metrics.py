"""Tests for molecular dynamics trajectory metrics."""

import itertools
import re
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from ase import Atoms, units
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.singlepoint import SinglePointCalculator

from matbench_discovery.metrics import md as md_metrics

np_rng = np.random.default_rng(seed=0)


class OffsetCalculator(Calculator):
    """Returns the frame's stored results shifted by constant offsets."""

    implemented_properties = ("energy", "forces")

    def __init__(
        self, energy_offset: float = 0, force_offset: float = 0, **kwargs: object
    ) -> None:
        """Store constant energy and force offsets applied to reference labels."""
        super().__init__(**kwargs)
        self.energy_offset = energy_offset
        self.force_offset = force_offset

    def calculate(
        self,
        atoms: Atoms | None = None,
        properties: list[str] | None = None,
        system_changes: list[str] = all_changes,
    ) -> None:
        """Add offsets to the reference energy/forces stored on the frame."""
        super().calculate(atoms, properties, system_changes)
        assert atoms is not None
        self.results["energy"] = atoms.info["ref_energy"] + self.energy_offset
        self.results["forces"] = atoms.arrays["ref_forces"] + self.force_offset


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
    """Store reference energy/forces where OffsetCalculator and ASE getters find
    them: in info/arrays and as SinglePointCalculator results.
    """
    atoms.info["ref_energy"] = energy
    atoms.arrays["ref_forces"] = forces
    atoms.calc = SinglePointCalculator(atoms, energy=energy, forces=forces)
    return atoms


def make_jiggled_frames(
    n_frames: int, *, pressure_gpa: float | None = None
) -> list[Atoms]:
    """Periodic 8-atom cubic-lattice frames with small random displacements."""
    base_positions = 3.0 * np.array(
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


@pytest.mark.parametrize(
    ("formulas", "n_bins", "err_msg"),
    [
        ([], 50, "empty trajectory"),
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


# === velocities and VDOS ===


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


def test_calc_pressure_voigt_vs_matrix() -> None:
    """Voigt and 3x3 stress representations should give identical pressures."""
    pressure_gpa = 2.5
    stress_matrix = -pressure_gpa * units.GPa * np.eye(3)
    stress_voigt = np.array([*np.diag(stress_matrix), 0, 0, 0])

    assert md_metrics.calc_pressure(stress_matrix) == pytest.approx(pressure_gpa)
    assert md_metrics.calc_pressure(stress_voigt) == pytest.approx(pressure_gpa)
    with pytest.raises(ValueError, match=re.escape("shape (6,) or (3, 3)")):
        md_metrics.calc_pressure(np.zeros(5))


def test_get_trajectory_pressures() -> None:
    """Per-frame pressures should be read from attached stress results."""
    pressures = [0.5, -1.0, 2.0]
    frames = [
        make_jiggled_frames(1, pressure_gpa=pressure)[0] for pressure in pressures
    ]
    np.testing.assert_allclose(
        md_metrics.get_trajectory_pressures(frames), pressures, rtol=1e-10, atol=1e-12
    )


def test_calc_pressure_metrics() -> None:
    """Offset distributions should give the offset; W1 must be order-independent."""
    p_ref = np.array([0.0, 1.0, 2.0, 3.0])
    p_shifted = p_ref + 0.5

    assert md_metrics.calc_pressure_metrics(p_ref, p_ref) == {
        "pressure_mae": 0,
        "pressure_wasserstein": 0,
    }
    shifted = md_metrics.calc_pressure_metrics(p_ref, p_shifted)
    assert shifted["pressure_mae"] == pytest.approx(0.5)
    assert shifted["pressure_wasserstein"] == pytest.approx(0.5)
    # unequal lengths truncate the MAE to the shorter trajectory
    truncated = md_metrics.calc_pressure_metrics(p_ref, p_shifted[:2])
    assert truncated["pressure_mae"] == pytest.approx(0.5)
    # W1 is invariant to frame order, MAE is not
    reversed_order = md_metrics.calc_pressure_metrics(p_ref, p_shifted[::-1])
    assert reversed_order["pressure_wasserstein"] == pytest.approx(0.5)
    assert reversed_order["pressure_mae"] > 0.5

    with pytest.raises(ValueError, match="empty"):
        md_metrics.calc_pressure_metrics(np.array([]), p_ref)

    # MAE must pair identical timestamps when cadences differ: pred samples ref
    # at every 2nd timestamp plus a pure +10 offset, so aligned MAE is exactly 10
    # (index pairing would wrongly give 11.5)
    p_ref_fine = np.arange(7.0)  # t = 0..6 at dt=1
    p_pred_coarse = p_ref_fine[::2] + 10  # t = 0,2,4,6 at dt=2
    aligned = md_metrics.calc_pressure_metrics(
        p_ref_fine, p_pred_coarse, ref_time_step_fs=1, pred_time_step_fs=2
    )
    assert aligned["pressure_mae"] == pytest.approx(10)

    # Awkward cadence ratio: ref dt=0.3, pred dt=0.7 -> common grid every 2.1 fs.
    # The common timestamps are 0 and 2.1 fs (ref indices 0,7; pred indices 0,3).
    p_ref_awkward = np.arange(8.0)
    p_pred_awkward = np.array([100.0, -1.0, -2.0, 107.0])
    awkward = md_metrics.calc_pressure_metrics(
        p_ref_awkward,
        p_pred_awkward,
        ref_time_step_fs=0.3,
        pred_time_step_fs=0.7,
    )
    assert awkward["pressure_mae"] == pytest.approx(100)


# === combined error ===


@pytest.mark.parametrize(
    ("rdf_error", "vdos_error", "expected"),
    [
        (0, 0, 0),
        (10, 0, 10),  # all weight on the worse metric
        (0, 10, 10),
        (10, 10, 10),
        (20, 10, 20 / 30 * 20 + 10 / 30 * 10),  # error-weighted, > plain mean of 15
    ],
)
def test_calc_combined_error(
    rdf_error: float, vdos_error: float, expected: float
) -> None:
    """Combined error should weight the worse of RDF and VDOS errors more."""
    assert md_metrics.calc_combined_error(rdf_error, vdos_error) == pytest.approx(
        expected
    )

    with pytest.raises(ValueError, match="must be non-negative"):
        md_metrics.calc_combined_error(-1, 5)


# === energy/force RMSE ===


@pytest.mark.parametrize("force_offset", [0, 0.1, 0.25])
def test_calc_energy_force_rmse(force_offset: float) -> None:
    """Raw per-atom energy RMSE preserves constant offsets; forces map directly."""
    frames = [
        attach_ref_labels(
            Atoms("H2", positions=np_rng.random((2, 3)) * 2, cell=[5] * 3),
            energy=float(np_rng.normal()),
            forces=np_rng.normal(size=(2, 3)),
        )
        for _ in range(3)
    ]
    calc = OffsetCalculator(energy_offset=0.8, force_offset=force_offset)
    result = md_metrics.calc_energy_force_rmse(frames, calc)
    assert result["energy_rmse"] == pytest.approx(0.4, abs=1e-12)
    assert result["force_rmse"] == pytest.approx(abs(force_offset), abs=1e-12)

    with pytest.raises(ValueError, match="empty trajectory"):
        md_metrics.calc_energy_force_rmse([], calc)
    single_frame = md_metrics.calc_energy_force_rmse(frames[:1], calc)
    assert single_frame["energy_rmse"] == pytest.approx(0.4, abs=1e-12)


def test_calc_energy_force_rmse_includes_offsets_and_fluctuations() -> None:
    """Energy RMSE should follow paper Eq. 1 without subtracting mean offsets."""
    n_atoms = 2
    fluctuations = np.array([0.0, 0.3, -0.3, 0.6, -0.6])  # eV total, varying
    forces = np.zeros((n_atoms, 3))
    frames = []
    for fluct in fluctuations:
        atoms = Atoms("H2", positions=np_rng.random((n_atoms, 3)) * 2, cell=[5] * 3)
        base = float(np_rng.normal())
        atoms.calc = SinglePointCalculator(atoms, energy=base, forces=forces)  # ref
        atoms.info["ref_energy"] = base + fluct  # OffsetCalculator adds +1000
        atoms.arrays["ref_forces"] = forces
        frames.append(atoms)

    result = md_metrics.calc_energy_force_rmse(frames, OffsetCalculator(1000))
    expected = np.sqrt(np.mean(((fluctuations + 1000) / n_atoms) ** 2))
    assert result["energy_rmse"] == pytest.approx(expected)


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
        "vdos_error",
        "pressure_mae",
        "pressure_wasserstein",
    }
    assert 0 <= metrics["rdf_error"] <= 100
    assert 0 <= metrics["vdos_error"] <= 100
    assert metrics["pressure_mae"] == pytest.approx(0.2, abs=1e-10)
    assert metrics["pressure_wasserstein"] == pytest.approx(0.2, abs=1e-10)

    # frames without stress data should yield NaN pressure metrics, not errors
    no_stress = md_metrics.evaluate_md_system(
        make_jiggled_frames(8),
        make_jiggled_frames(8),
        ref_time_step_fs=1,
        pred_time_step_fs=1,
    )
    assert np.isnan(no_stress["pressure_mae"])
    assert np.isnan(no_stress["pressure_wasserstein"])


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
        calculator=OffsetCalculator(energy_offset=0.8),
    )
    assert metrics["energy_rmse"] == pytest.approx(0.1)  # 0.8 eV / 8 atoms
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
            "vdos_error": [30.0, 10.0],
            "pressure_mae": [1.0, np.nan],  # one system without stress data
            "pressure_wasserstein": [0.5, np.nan],
        }
    )
    metrics = md_metrics.calc_md_metrics(df_md)

    assert metrics["energy_rmse"] == pytest.approx(0.2)
    assert metrics["force_rmse"] == pytest.approx(0.3)
    assert metrics["rdf_error"] == pytest.approx(15)
    assert metrics["vdos_error"] == pytest.approx(20)
    assert metrics["pressure_mae"] == pytest.approx(1)
    assert metrics["pressure_wasserstein"] == pytest.approx(0.5)
    assert metrics["combined_error"] == pytest.approx(
        md_metrics.calc_combined_error(15, 20)
    )
    assert metrics["n_systems"] == 2

    # partial columns work, unrecognized columns raise
    partial = md_metrics.calc_md_metrics(pd.DataFrame({"rdf_error": [10.0, 20.0]}))
    assert partial == {"rdf_error": 15, "n_systems": 2}
    with pytest.raises(ValueError, match="No recognized MD metric columns"):
        md_metrics.calc_md_metrics(pd.DataFrame({"unrelated": [1]}))


def test_load_per_system_metrics(tmp_path: Path) -> None:
    """Per-system CSVs concat into one frame indexed by system; reruns override."""
    (tmp_path / "sysA.csv").write_text("system,rdf_error\nsysA,10.0\n")
    (tmp_path / "sysB.csv").write_text("system,rdf_error\nsysB,20.0\n")
    # a rerun of sysA with a better value, written later -> must override the first
    (tmp_path / "sysA_rerun.csv").write_text("system,rdf_error\nsysA,5.0\n")

    paths = [f"{tmp_path}/{name}.csv" for name in ("sysA", "sysB", "sysA_rerun")]
    df_md = md_metrics.load_per_system_metrics(paths)
    assert list(df_md.index) == ["sysB", "sysA"]  # sysA dedup'd to its last (rerun) row
    assert df_md.loc["sysA", "rdf_error"] == 5.0
    assert df_md.loc["sysB", "rdf_error"] == 20.0

    with pytest.raises(ValueError, match="No per-system metric CSVs given"):
        md_metrics.load_per_system_metrics([])
    (tmp_path / "no_sys.csv").write_text("rdf_error\n10.0\n")
    with pytest.raises(ValueError, match="lack a 'system' column"):
        md_metrics.load_per_system_metrics([f"{tmp_path}/no_sys.csv"])


@pytest.mark.parametrize(
    "init_yaml",  # placeholder strings at leaf and intermediate level get replaced
    ["metrics:\n  md: not available\n", "metrics: not available\n"],
)
def test_write_metrics_to_yaml(tmp_path: Path, init_yaml: str) -> None:
    """Metrics should be written under metrics.md with units and rounding."""
    yaml_path = tmp_path / "model.yml"
    yaml_path.write_text(f"model_name: test\n{init_yaml}")
    model = SimpleNamespace(yaml_path=str(yaml_path))

    metrics = {
        "energy_rmse": 0.123456,
        "force_rmse": 0.234567,
        "rdf_error": 12.34567,
        "combined_error": 23.456789,
        "n_systems": 17,
    }
    md_metrics.write_metrics_to_yaml(
        model,  # ty: ignore[invalid-argument-type]
        metrics,
        pred_file_path="models/test/md-metrics.csv",
        pred_file_url="https://example.com/md-metrics.csv",
    )

    text = yaml_path.read_text()
    assert "pred_file: models/test/md-metrics.csv" in text
    assert "pred_file_url: https://example.com/md-metrics.csv" in text
    assert "energy_rmse: 0.1235 # eV/atom" in text
    assert "force_rmse: 0.2346 # eV/Å" in text
    assert "rdf_error: 12.3457 # %" in text
    assert "combined_error: 23.4568 # %" in text
    assert "n_systems: 17 # count" in text
    # the 'not available' placeholder must be replaced, not kept
    assert "not available" not in text
