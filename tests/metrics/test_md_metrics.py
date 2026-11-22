"""Tests for molecular dynamics trajectory metrics."""

import itertools
import re
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml
from matscikit import FEMTOSECOND, K_B, Structure
from matscikit.elements import ATOMIC_MASSES
from matscikit.io import Frame

from matbench_discovery.data import update_yaml_file
from matbench_discovery.enums import Model
from matbench_discovery.metrics import md as md_metrics
from tests.utils import patch_model_metadata

ANISOTROPIC_STRESS = -np.diag([1.0, 2.0, 6.0]) / md_metrics.EV_PER_A3_TO_GPA


@pytest.fixture
def np_rng() -> np.random.Generator:
    """Fresh seeded generator per test so inputs don't depend on test order."""
    return np.random.default_rng(seed=0)


def make_frame(
    formula: str,
    positions: object = None,
    *,
    cell: object,
    pbc: bool = False,
    **results: object,
) -> Frame:
    """ASE-Atoms-like frame constructor: formula like "CH4", cubic edge list or 3x3
    cell, scalar pbc, and optional energy/forces/stress results.
    """
    symbols = [
        symbol
        for symbol, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula)
        for _ in range(int(count or 1))
    ]
    cell_array = np.asarray(cell, dtype=float)
    lattice = np.diag(cell_array) if cell_array.ndim == 1 else cell_array
    coords = np.zeros((len(symbols), 3)) if positions is None else positions
    structure = Structure(lattice, symbols, coords, pbc=(pbc,) * 3)
    return Frame(structure, **results)  # ty: ignore[invalid-argument-type]


class ConstantModel:
    """Predicts a fixed energy and uniform force regardless of geometry."""

    def __init__(self, energy: float = 0, force: float = 0) -> None:
        """Store the constant energy and per-component force to return."""
        self.energy = energy
        self.force = force

    def compute(
        self, structures: Sequence[Structure], properties: Sequence[str]
    ) -> list[dict[str, object]]:
        """Return the stored constant energy and a uniform force on every atom."""
        del properties
        return [
            {"energy": self.energy, "forces": np.full((len(struct), 3), self.force)}
            for struct in structures
        ]


def make_jiggled_frames(
    n_frames: int, np_rng: np.random.Generator, *, pressure_gpa: float | None = None
) -> list[Frame]:
    """Periodic 8-atom cubic-lattice frames with small random displacements. The 2.4 A
    spacing keeps Si-Si pairs within the covalent-radius bond cutoff so matscikit.adf
    finds bonded triplets (real Si-Si bond ~2.35 A; covalent cutoff ~2.78 A).
    """
    base_positions = 2.4 * np.array(
        list(itertools.product(range(2), repeat=3)), dtype=float
    )
    frames = []
    for _ in range(n_frames):
        positions = base_positions + np_rng.normal(scale=0.05, size=(8, 3))
        stress = (
            None
            if pressure_gpa is None
            else -pressure_gpa / md_metrics.EV_PER_A3_TO_GPA * np.eye(3)
        )
        frames.append(
            make_frame("Si8", positions, cell=[6] * 3, pbc=True, stress=stress)
        )
    return frames


def h2_frame(
    energy: float,
    np_rng: np.random.Generator,
    *,
    positions: np.ndarray | None = None,
    forces: np.ndarray | None = None,
) -> Frame:
    """An H2 reference frame carrying private energy/force labels."""
    if positions is None:
        positions = np_rng.random((2, 3)) * 2
    if forces is None:
        forces = np.zeros((2, 3))
    return make_frame("H2", positions, cell=[5] * 3, energy=energy, forces=forces)


# === RDF ===


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


@pytest.mark.parametrize(
    ("g_r_ref", "g_r_pred", "err_msg"),
    [
        (np.ones(3), np.ones(4), "differ"),  # mismatched grid lengths
        (np.array([0.0, -1.0, 1.0]), np.ones(3), "negative"),
        (np.array([0.0, np.inf, 1.0]), np.ones(3), "non-finite"),
    ],
)
def test_calc_rdf_error_invalid_input(
    g_r_ref: np.ndarray, g_r_pred: np.ndarray, err_msg: str
) -> None:
    """RDF error should reject mismatched, negative and non-finite distributions."""
    with pytest.raises(ValueError, match=err_msg):
        md_metrics.calc_rdf_error(np.array([1.0, 2.0, 3.0]), g_r_ref, g_r_pred)


# === ADF ===


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


# === trajectory observables, velocities and vDOS ===


def test_trajectory_observables(np_rng: np.random.Generator) -> None:
    """Frame lists coerce to the matscikit kernels' arrays on the default r_max grid."""
    observables = md_metrics.trajectory_observables(
        make_jiggled_frames(8, np_rng), time_step_fs=2, n_rdf_bins=50
    )
    assert {key: len(grid) for key, (grid, _) in observables.items()} == {
        "rdf": 50,
        "adf": 180,
        "vdos": 5,
    }
    # 8 atoms in a 6 A cube: r_max = 3 A, Si-Si bonds within 2.78 A
    assert observables["rdf"][0][-1] == pytest.approx(3 - 3 / 100, rel=1e-14)
    np.testing.assert_allclose(observables["adf"][1].sum(), 1, rtol=1e-12)


def test_equipartition_temperature() -> None:
    """Equipartition T recovers the velocities' temperature and scales as 1/dt**2 (the
    property that catches a mislabeled dt_fs: e.g. dt 5x too large -> T 25x too low).
    """
    rng = np.random.default_rng(0)
    n_atoms, n_frames, dt = 32, 6, 1.5  # fs
    velocity = rng.normal(scale=0.01, size=(n_atoms, 3))  # A/fs, constant per atom
    base = rng.uniform(0, 10, size=(n_atoms, 3))
    frames = [
        make_frame(
            "Cu32", positions=base + velocity * (step * dt), cell=[80] * 3, pbc=False
        )
        for step in range(n_frames)
    ]
    temp = md_metrics.equipartition_temperature(frames, time_step_fs=dt)
    # independent expectation from the same constant velocities (linear motion -> exact)
    # 1 amu (A/fs)^2 in eV in matscikit's units (its CODATA, not a hardcoded older one)
    ke_ev = 0.5 * ATOMIC_MASSES[29] * (velocity**2).sum() * FEMTOSECOND**-2
    expected = 2 * ke_ev / (3 * n_atoms * K_B)
    assert temp == pytest.approx(expected, rel=1e-9)
    # halving the assumed dt quadruples the inferred temperature
    half_dt = md_metrics.equipartition_temperature(frames, time_step_fs=dt / 2)
    assert half_dt == pytest.approx(4 * temp, rel=1e-9)


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
    ("vdos_ref", "band_frac", "err_msg"),
    [
        (np.full(11, -1.0), 0.5, "negative intensities"),
        (np.zeros(11), 0.5, "non-positive area"),
        (np.full(11, np.inf), 0.5, "non-finite intensities"),
        *((np.ones(11), band_frac, "band_frac") for band_frac in (0.0, -0.1, 1.5)),
    ],
)
def test_calc_vdos_error_invalid_input(
    vdos_ref: np.ndarray, band_frac: float, err_msg: str
) -> None:
    """Reject negative, zero-area or non-finite spectra and band_frac outside (0, 1]."""
    freqs = np.linspace(0, 10, 11)
    with pytest.raises(ValueError, match=err_msg):
        md_metrics.calc_vdos_error(
            freqs, vdos_ref, freqs, np.ones(11), band_frac=band_frac
        )


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


def test_get_trajectory_pressures(np_rng: np.random.Generator) -> None:
    """Pressures average the full 3D stress trace, not a 2D slab formula, so the
    anisotropic frame gives 3 GPa, not 1.5.
    """
    pressures = [0.5, -1.0, 2.0]
    frames = [
        make_jiggled_frames(1, np_rng, pressure_gpa=pressure)[0]
        for pressure in pressures
    ]
    aniso = make_frame(
        "Si8",
        np_rng.normal(scale=0.05, size=(8, 3)),
        cell=[6] * 3,
        pbc=True,
        stress=ANISOTROPIC_STRESS,
    )
    np.testing.assert_allclose(
        md_metrics.get_trajectory_pressures([*frames, aniso]),
        [*pressures, 3.0],
        rtol=1e-10,
        atol=1e-12,
    )


@pytest.mark.parametrize(
    ("p_ref", "p_pred", "expected_mae", "expected_w1"),
    [
        (np.arange(4.0), np.arange(4.0), 0.0, 0.0),
        (np.arange(4.0), np.arange(4.0) + 0.5, 0.5, 0.5),
        # Unequal lengths: mean bias 0.5; W1 = (0.5 + 0.5 + 0.5 + 1.5) / 4.
        (np.arange(4.0), np.array([0.5, 1.5]), 0.5, 0.75),
        (np.arange(7.0), np.arange(7.0)[::-1], 0.0, 0.0),
        (np.arange(4.0), (np.arange(4.0) + 0.5)[::-1], 0.5, 0.5),
    ],
    ids=["identical", "offset", "unequal-lengths", "reversed", "reversed-offset"],
)
def test_calc_pressure_metrics(
    p_ref: np.ndarray, p_pred: np.ndarray, expected_mae: float, expected_w1: float
) -> None:
    """Pressure bias and Wasserstein distance allow reordered and unequal series."""
    metrics = md_metrics.calc_pressure_metrics(p_ref, p_pred)
    assert metrics["pressure_mae"] == pytest.approx(expected_mae)
    assert metrics["pressure_wasserstein"] == pytest.approx(
        expected_w1, rel=1e-6, abs=0
    )


# === private energy/force diagnostics ===


def test_calc_energy_force_rmse_uses_private_labels(
    np_rng: np.random.Generator,
) -> None:
    """Energy fluctuations and force RMSE come from private reference labels."""
    frames = [
        h2_frame(1.0, np_rng, forces=np.zeros((2, 3))),
        h2_frame(3.0, np_rng, forces=np.ones((2, 3))),
    ]

    metrics = md_metrics.calc_energy_force_rmse(
        frames, ConstantModel(energy=2.0, force=0.0)
    )

    # predicted energy fluctuations are zero; reference fluctuations are ±0.5 eV/atom
    assert metrics["energy_rmse"] == pytest.approx(0.5)
    assert metrics["force_rmse"] == pytest.approx(np.sqrt(0.5))


@pytest.mark.parametrize(
    ("results", "message"), [(1, "2 model results"), (2, "forces of shape")]
)
def test_energy_force_rmse_rejects_short_batch(
    results: int,
    message: str,
    np_rng: np.random.Generator,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Missing frames or atom forces cannot broadcast into a valid RMSE."""
    frames = [h2_frame(1.0, np_rng)] * 2
    model = ConstantModel()
    monkeypatch.setattr(
        model,
        "compute",
        lambda *_: [{"energy": 1.0, "forces": np.zeros((1, 3))}] * results,
    )
    with pytest.raises(ValueError, match=message):
        md_metrics.calc_energy_force_rmse(frames, model)


def test_energy_force_rmse_rejects_missing_private_labels(
    np_rng: np.random.Generator,
) -> None:
    """Private diagnostics must fail loudly when the reference lacks labels."""
    atoms = make_frame("H2", positions=np_rng.random((2, 3)) * 2, cell=[5] * 3)

    with pytest.raises(ValueError, match="energy/forces"):
        md_metrics.calc_energy_force_rmse([atoms], ConstantModel())


@pytest.mark.parametrize(
    ("shift", "expected"),
    [(0, 0), (1000, 100)],  # identical -> 0% non-overlap; far apart -> ~100%
    ids=["identical", "disjoint"],
)
def test_calc_pressure_histogram_error(
    shift: float, expected: float, np_rng: np.random.Generator
) -> None:
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


@pytest.mark.parametrize(
    ("func_name", "p_ref", "kwargs", "err_msg"),
    [
        ("calc_pressure_metrics", np.array([]), {}, "empty"),
        ("calc_pressure_histogram_error", np.array([]), {}, "empty"),
        ("calc_pressure_histogram_error", np.arange(4.0), {"n_bins": 1}, "bins"),
        ("calc_pressure_histogram_error", np.array([0.0, np.nan]), {}, "finite"),
    ],
)
def test_pressure_metrics_reject_bad_input(
    func_name: str, p_ref: np.ndarray, kwargs: dict[str, int], err_msg: str
) -> None:
    """Empty arrays, <2 bins and non-finite pressures raise, not return garbage."""
    with pytest.raises(ValueError, match=err_msg):
        getattr(md_metrics, func_name)(p_ref, np.arange(4.0), **kwargs)


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


@pytest.mark.parametrize(("ref_dt", "pred_dt"), [(0, 1), (1, 0)])
def test_matched_frame_counts_rejects_zero_time_step(
    ref_dt: float, pred_dt: float
) -> None:
    """A zero time step has no common grid, so it must raise rather than divide by 0."""
    with pytest.raises(ValueError, match="must be positive"):
        md_metrics.matched_frame_counts(
            n_ref_frames=5,
            n_pred_frames=5,
            ref_time_step_fs=ref_dt,
            pred_time_step_fs=pred_dt,
        )


def test_evaluate_md_system_rejects_short_trajectories(
    np_rng: np.random.Generator,
) -> None:
    """Trajectories shorter than two matched frames should raise."""
    with pytest.raises(ValueError, match="too short after time matching"):
        md_metrics.evaluate_md_system(
            make_jiggled_frames(1, np_rng),
            make_jiggled_frames(10, np_rng),
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
            "load_time_sec": [3.0, 0.0],  # one process ran both systems
            "n_processes": [1, 0],
            "max_rss_gb": [4.2, 6.1],
            "max_gpu_mem_gb": [11.5, 9.0],
            "hardware": ["NVIDIA H200", "NVIDIA H200"],
            "cost_source": ["runner"] * 2,
            # md.run_info_row CSV encoding of the acceleration list
            "acceleration": ["cuequivariance"] * 2,
        }
    )
    metrics = md_metrics.calc_md_metrics(df_md)

    labels = {key: metrics.pop(key) for key in ("cost_source", "hardware")}
    assert labels == {"cost_source": "runner", "hardware": "NVIDIA H200"}
    assert metrics.pop("acceleration") == ["cuequivariance"]
    # exact-key equality also proves combined_score is absent (CMDS is site-computed)
    assert metrics == pytest.approx(
        {
            "run_time_sec": 300.75,  # summed, not averaged
            "load_time_sec": 3.0,
            "n_processes": 1,
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
    # "none" decodes to no acceleration
    no_accel = df_md.assign(acceleration="none")
    assert md_metrics.calc_md_metrics(no_accel)["acceleration"] == []
    # a partial label (e.g. an untagged cached pre-provenance rollout) drops the
    # whole cost block, never the metrics
    with pytest.warns(UserWarning, match="mixed or partial cost_source"):
        mixed = md_metrics.calc_md_metrics(df_md.assign(cost_source=["runner", None]))
    assert set(mixed).isdisjoint(md_metrics.MD_COST_COLS)
    assert mixed["rdf_error"] == pytest.approx(15)


@pytest.mark.parametrize(
    ("cols", "dropped"),
    [
        (
            {"energy_rmse": [0.001, np.nan], "force_rmse": [0.1, np.nan]},
            ("energy_rmse", "force_rmse"),
        ),
        (
            # model without stress support: no pressure mean to publish (NaN would
            # land in the YAML), while non-pressure metrics stay
            {
                "pressure_mae": [np.nan, np.nan],
                "pressure_wasserstein": [np.nan, np.nan],
                "pressure_error": [np.nan, np.nan],
            },
            ("pressure_mae", "pressure_wasserstein", "pressure_error"),
        ),
        (
            # pressure columns disagreeing on which systems are NaN -> inconsistent
            # per-system coverage, drop all pressure metrics
            {
                "pressure_mae": [1.0, 2.0],
                "pressure_wasserstein": [0.5, np.nan],
                "pressure_error": [40.0, 30.0],
            },
            ("pressure_mae", "pressure_wasserstein", "pressure_error"),
        ),
    ],
    ids=[
        "partial_energy_force",
        "all_nan_pressure",
        "mismatched_pressure_nans",
    ],
)
def test_calc_md_metrics_all_or_nothing_columns(
    cols: dict[str, list[float]], dropped: tuple[str, ...]
) -> None:
    """E/F and pressure means publish only with full per-system coverage."""
    df_md = pd.DataFrame({"rdf_error": [10.0, 20.0]} | cols)
    metrics = md_metrics.calc_md_metrics(df_md)
    assert set(dropped).isdisjoint(metrics)
    assert set(cols) - set(dropped) <= set(metrics)
    assert metrics["rdf_error"] == pytest.approx(15)
    assert metrics["n_systems"] == 2


def test_write_metrics_to_yaml(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """metrics.md gets units, rounding and the artifact's cost after pred_file."""
    yaml_path = tmp_path / "model.yml"
    yaml_path.write_text("model_name: test\nmetrics: {}\n", encoding="utf-8")
    model = Model["mace_mp_0"]
    patch_model_metadata(monkeypatch, model, yaml_path=str(yaml_path))

    cost = {
        "cost_source": "runner",
        "hardware": "NVIDIA H200",
        "acceleration": ["cuequivariance"],
        "run_time_sec": 4521.484,
        "load_time_sec": 12.5,
    }
    metrics: dict[str, object] = {
        **cost,
        "energy_rmse": 1.23456,
        "force_rmse": 98.7654,
        "rdf_error": 12.34567,
        "n_systems": 17,
    }
    path = "models/mace/mace-mp-0/2026-07-01-md-metrics.csv.gz"
    url = "https://example.com/2026-07-01-md-metrics.csv.gz"
    md_metrics.write_metrics_to_yaml(
        model, metrics, pred_file_path=path, pred_file_url=url
    )

    text = yaml_path.read_text(encoding="utf-8")
    md = yaml.safe_load(text)["metrics"]["md"]
    assert list(md)[:2] == ["pred_file", "cost"]
    assert md["pred_file"] == {"name": path, "url": url}
    assert "  cost: {source: runner, hardware: NVIDIA H200, acceleration: " in text
    assert "energy_rmse: 1.2346 # meV/atom" in text
    assert "force_rmse: 98.7654 # meV/Å" in text
    assert "rdf_error: 12.3457 # %" in text
    assert "n_systems: 17 # count" in text
    # CMDS is site-computed, never persisted
    assert "combined_score" not in text

    # recomputing the registered CSV re-derives the whole cost block from it: a
    # recompute whose CSV has no cost contradicts the runner tag and raises
    recomputed = md_metrics.write_metrics_to_yaml(
        model, {**cost, "run_time_sec": 10.0, "rdf_error": 11.0, "n_systems": 17}
    )
    assert recomputed["cost"]["run_time_sec"] == 10.0
    assert recomputed["pred_file"] == {"name": path, "url": url}
    with pytest.raises(ValueError, match="records no cost"):
        md_metrics.write_metrics_to_yaml(model, {"rdf_error": 11.0, "n_systems": 17})

    # a new source without cost columns cannot inherit the old run's cost
    md_metrics.write_metrics_to_yaml(
        model, {"rdf_error": 3.0, "n_systems": 17}, pred_file_path=path
    )
    assert yaml.safe_load(yaml_path.read_text())["metrics"]["md"] == {
        "pred_file": {"name": path},
        "rdf_error": 3.0,
        "n_systems": 17,
    }

    # an untagged legacy CSV cost with no YAML cost block is skipped, not published
    legacy_cost = {"hardware": "NVIDIA H100", "run_time_sec": 7.0, "max_rss_gb": 2.4}
    recomputed = md_metrics.write_metrics_to_yaml(
        model, {**legacy_cost, "rdf_error": 3.0, "n_systems": 17}
    )
    assert "cost" not in recomputed
    assert "skipping MD cost" in capsys.readouterr().out

    # a recompute keeps the reproducibility verdict, registering a new CSV drops it
    repro = {"status": "reproduced", "rel_changes": {"rdf_error": 0.0}}
    update_yaml_file(str(yaml_path), "metrics.md", {"reproducibility": repro})
    scores = {"rdf_error": 3.0, "n_systems": 17}
    assert md_metrics.write_metrics_to_yaml(model, scores)["reproducibility"] == repro
    written = md_metrics.write_metrics_to_yaml(model, scores, pred_file_path=path)
    assert "reproducibility" not in written
