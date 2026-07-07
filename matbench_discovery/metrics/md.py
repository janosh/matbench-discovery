"""Metrics comparing MLIP molecular dynamics trajectories to ab-initio references.

Implements the finite-temperature observables proposed for the MD benchmark task
(https://github.com/janosh/matbench-discovery/pull/344): radial distribution
function (RDF) error, angular distribution function (ADF) error, pressure metrics
from the stress tensor trace, vibrational density of states (vDOS) error from the
velocity autocorrelation function, and private-label energy/force diagnostics.
The combined MD score (CMDS) is computed on the fly by the site from these
components (like CPS), never here or in model YAMLs.
"""

import math
from collections.abc import Sequence
from fractions import Fraction
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from ase import Atoms, units
from ase.data import atomic_masses, covalent_radii
from ase.geometry import find_mic, get_distances, minkowski_reduce
from ruamel.yaml.comments import CommentedMap

from matbench_discovery import repo_relative_path
from matbench_discovery.data import update_yaml_file
from matbench_discovery.trajectory import Trajectory

if TYPE_CHECKING:
    from ase.calculators.calculator import Calculator

    from matbench_discovery.enums import Model

TrajectoryLike = Trajectory | Sequence[Atoms]


def _as_trajectory(trajectory: TrajectoryLike) -> Trajectory:
    """Coerce a Trajectory or a sequence of ASE Atoms to a Trajectory, so metrics run
    on dense arrays regardless of input. Passing a Trajectory (the hot path) is a no-op.
    """
    if isinstance(trajectory, Trajectory):
        return trajectory
    return Trajectory.from_ase(list(trajectory))


def _check_density(values: np.ndarray, label: str, *, check_area: bool = True) -> None:
    """Reject non-finite, negative or (when check_area) zero-area distribution densities
    shared by the RDF/ADF/vDOS error metrics.
    """
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{label} has non-finite intensities")
    if np.any(values < 0):
        raise ValueError(f"{label} has negative intensities")
    if check_area and values.sum() <= 0:
        raise ValueError(f"{label} has non-positive area")


# units of per-system metric columns shared between eval scripts and aggregation
# (no combined_score here: CMDS is site-computed, see module docstring)
# NOTE: these are the *model-level* units after calc_md_metrics aggregation. In the
# per-system CSVs written by run_md_benchmark, energy_rmse/force_rmse are still in
# eV/atom and eV/Å (calc_energy_force_rmse output); calc_md_metrics scales them x1000.
METRIC_UNITS = {
    "run_time_sec": "s",
    "max_rss_gb": "GB",
    "max_gpu_mem_gb": "GB",
    "energy_rmse": "meV/atom",
    "force_rmse": "meV/Å",
    "rdf_error": "%",
    "adf_error": "%",
    "vdos_error": "%",
    "pressure_mae": "GPa",
    "pressure_wasserstein": "GPa",
    "pressure_error": "%",  # pressure-histogram overlap error (paper Eq. 9)
    "n_systems": "count",
}
# run-provenance columns aggregated by sum/max/passthrough instead of the mean
RUN_PROVENANCE_COLS = ("run_time_sec", "max_rss_gb", "max_gpu_mem_gb", "hardware")
# mean-aggregated per-system error columns
PER_SYSTEM_METRIC_COLS = tuple(
    key for key in METRIC_UNITS if key != "n_systems" and key not in RUN_PROVENANCE_COLS
)


def min_image_radius(cells: np.ndarray, pbc: np.ndarray) -> float:
    """Largest radius at which minimum-image pair distances count each pair at most
    once = half the shortest lattice vector (Minkowski first successive minimum),
    minimized over frames.

    Minimum-image distances record only the nearest periodic image of each atom, so a
    histogram (RDF) is complete only while no two images of an atom lie within the
    cutoff sphere; two distinct images are at least one shortest lattice vector apart,
    giving the half-shortest-vector limit. For skewed or non-reduced cells a lattice-
    vector combination can be shorter than every basis vector, so half the shortest
    *basis* vector overestimates this and silently corrupts g(r) at large r.
    """

    def shortest_lattice_vector(cell: np.ndarray) -> float:
        reduced, _ = minkowski_reduce(cell, pbc)
        lengths = np.linalg.norm(reduced, axis=1)
        periodic = lengths[pbc]  # only periodic directions are true lattice vectors
        return float(periodic.min() if periodic.size else lengths.min())

    # cells are usually constant (NVT) -> reduce once; else minimize over frames
    unique_cells = cells[:1] if np.allclose(cells, cells[0]) else cells
    return min(shortest_lattice_vector(cell) for cell in unique_cells) / 2


def calc_rdf(
    trajectory: TrajectoryLike,
    *,
    n_bins: int = 500,
    r_max: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Time-averaged all-pair radial distribution function g(r) of a trajectory.

    Returns bin-center distances r in Angstrom and g(r), histogrammed into n_bins
    between 0 and r_max. r_max defaults to half the smallest cell length of the
    first frame, the validity limit of minimum-image distances.
    """
    traj = _as_trajectory(trajectory)
    if traj.n_frames == 0:
        raise ValueError("Cannot compute RDF of empty trajectory")
    if n_bins < 2:
        raise ValueError(f"Need >= 2 RDF bins, got {n_bins=}")
    n_atoms = traj.n_atoms
    if n_atoms < 2:
        raise ValueError(f"Cannot compute RDF with {n_atoms=} < 2")
    mic_radius = min_image_radius(traj.cell, traj.pbc)
    if r_max is None:
        r_max = mic_radius
    elif r_max > mic_radius + 1e-12:
        raise ValueError(
            f"{r_max=} exceeds minimum-image validity limit {mic_radius:.6g} A "
            "for at least one trajectory frame"
        )
    if r_max <= 0:
        raise ValueError(f"{r_max=} must be positive")

    bin_edges = np.linspace(0, r_max, n_bins + 1)
    counts = np.zeros(n_bins)
    upper_tri_idx = np.triu_indices(n_atoms, k=1)
    inv_volume_sum = 0.0  # per-frame densities make this exact for varying cells
    for idx in range(traj.n_frames):
        cell = traj.cell[idx]
        # get_distances is the exact kernel Atoms.get_all_distances(mic=True) uses, so
        # results match the previous Atoms-based path to machine precision
        _, dist_matrix = get_distances(traj.positions[idx], cell=cell, pbc=traj.pbc)
        counts += np.histogram(dist_matrix[upper_tri_idx], bins=bin_edges)[0]
        inv_volume_sum += 1 / abs(float(np.linalg.det(cell)))

    radii = (bin_edges[:-1] + bin_edges[1:]) / 2
    shell_volumes = 4 / 3 * np.pi * np.diff(bin_edges**3)
    # Paper Eq. 4 normalization: rho*N shell volume. The upper-triangle histogram
    # counts unordered pairs, hence the /2 relative to Eq. 4's ordered i != j sum.
    ideal_counts = n_atoms**2 / 2 * inv_volume_sum * shell_volumes
    return radii, counts / ideal_counts


def calc_rdf_error(
    radii: np.ndarray, g_r_ref: np.ndarray, g_r_pred: np.ndarray
) -> float:
    """RDF error in percent (paper Eq. 3) between two g(r) on a shared radii grid.

        error = 100 * min(1, integral |g_ref - g_pred| dr / integral |g_ref - 1| dr)

    0 means perfect match, 100 means as different from the reference as an ideal
    gas (also returned if the reference itself is indistinguishable from an ideal
    gas, i.e. zero denominator).
    """
    if not len(radii) == len(g_r_ref) == len(g_r_pred):
        raise ValueError(f"{len(radii)=}, {len(g_r_ref)=}, {len(g_r_pred)=} differ")
    g_r_ref = np.asarray(g_r_ref, dtype=float)
    g_r_pred = np.asarray(g_r_pred, dtype=float)
    _check_density(g_r_ref, "reference RDF", check_area=False)
    _check_density(g_r_pred, "predicted RDF", check_area=False)
    denominator = np.trapezoid(np.abs(g_r_ref - 1), x=radii)
    if denominator == 0:
        return 100.0
    numerator = np.trapezoid(np.abs(g_r_ref - g_r_pred), x=radii)
    return float(min(1, numerator / denominator) * 100)


def calc_adf(
    trajectory: TrajectoryLike,
    *,
    n_bins: int = 180,
    bond_tolerance: float = 1.25,
) -> tuple[np.ndarray, np.ndarray]:
    """Time-averaged bond-angle distribution function of a trajectory.

    Returns bin-center angles in degrees and a unit-area angular density over
    0-180 degrees. For each frame and central atom, every unordered pair of *bonded*
    neighbors contributes one angle. A pair (i, j) counts as bonded when its minimum-
    image distance is within ``bond_tolerance * (r_cov[Z_i] + r_cov[Z_j])`` (ASE
    covalent radii). Species-aware cutoffs capture first-coordination-shell bonding
    angles (tetrahedral, octahedral, aromatic, ...) across mixed bond length scales
    (e.g. Pb-Br vs C-H in a hybrid perovskite). A single global cutoff instead lets the
    many long-range neighbor pairs dominate, collapsing every structure's ADF onto the
    universal sin(theta) background and washing out the bonding-angle signal.
    """
    traj = _as_trajectory(trajectory)
    if traj.n_frames == 0:
        raise ValueError("Cannot compute ADF of empty trajectory")
    if n_bins < 2:
        raise ValueError(f"Need >= 2 ADF bins, got {n_bins=}")
    if bond_tolerance <= 0:
        raise ValueError(f"{bond_tolerance=} must be positive")
    n_atoms = traj.n_atoms
    if n_atoms < 3:
        raise ValueError(f"Cannot compute ADF with {n_atoms=} < 3")
    # Clamp only periodic systems to the minimum-image limit; isolated molecules have
    # no image ambiguity, so keep their full covalent-radius cutoff.
    r_cov = covalent_radii[traj.atomic_numbers]
    bond_cutoffs = bond_tolerance * (r_cov[:, None] + r_cov[None, :])
    if traj.pbc.any():
        bond_cutoffs = np.minimum(bond_cutoffs, min_image_radius(traj.cell, traj.pbc))

    bin_edges = np.linspace(0, 180, n_bins + 1)
    counts = np.zeros(n_bins)
    for frame_idx in range(traj.n_frames):
        mic_vectors, distances = get_distances(
            traj.positions[frame_idx], cell=traj.cell[frame_idx], pbc=traj.pbc
        )
        bonded = (distances > 1e-12) & (distances <= bond_cutoffs)
        # directed bonds (center -> neighbor); np.nonzero is row-major so ``centers``
        # is non-decreasing and each center's neighbors form one contiguous run
        centers, neighbors = np.nonzero(bonded)
        bond_unit_vecs = (
            mic_vectors[centers, neighbors] / distances[centers, neighbors][:, None]
        )
        # all unordered neighbor pairs per center, computed only over bonded pairs (so
        # this is O(n * coord^2), not the O(n^3) of a full all-triplet cosine tensor)
        seg_starts = np.unique(centers, return_index=True)[1]
        seg_counts = np.diff(np.append(seg_starts, len(centers)))
        left_idx, right_idx = [], []
        for start, count in zip(seg_starts, seg_counts, strict=True):
            if count >= 2:  # need >= 2 neighbors to form an angle
                tri_a, tri_b = np.triu_indices(count, k=1)
                left_idx.append(start + tri_a)
                right_idx.append(start + tri_b)
        if not left_idx:
            continue
        left = np.concatenate(left_idx)
        right = np.concatenate(right_idx)
        cosines = np.einsum("pd,pd->p", bond_unit_vecs[left], bond_unit_vecs[right])
        pair_angles = np.degrees(np.arccos(np.clip(cosines, -1, 1)))
        counts += np.histogram(pair_angles, bins=bin_edges)[0]

    if not counts.any():
        raise ValueError("Cannot compute ADF: no bonded neighbor angle pairs found")

    angles = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_widths = np.diff(bin_edges)
    return angles, counts / np.sum(counts * bin_widths)


def calc_adf_error(
    angles: np.ndarray, adf_ref: np.ndarray, adf_pred: np.ndarray
) -> float:
    """Bond-angle distribution error in percent, the ADF analog of the RDF metric. A
    Matbench addition (the paper's structural/dynamical observables are RDF, pressure
    and vDOS only).

        error = 100 * min(1, W1(adf_ref, adf_pred) / W1(adf_ref, sin-background))

    W1 is the Wasserstein-1 (earth-mover) distance in degrees between the two area-
    normalized angular densities. Normalizing by the reference's own distance to the
    structureless ``sin(theta)`` random-neighbor background (the angular 'ideal gas',
    the angular analog of g(r) = 1) makes this a fractional structural error:
    0 = perfect match, 100 = as different from the reference as a featureless
    distribution (also returned when the reference is itself indistinguishable from
    the background).

    Wasserstein distance, not the L1 overlap used for RDF/pressure, is deliberate: it
    grows with *how far* bonding-angle peaks are shifted (e.g. tetrahedral 109.5 deg ->
    square-planar 90 deg) instead of saturating once peaks stop overlapping. Its lack
    of a histogram-binning noise floor keeps resolving the small angular errors that
    separate otherwise-good MLIPs.
    """
    from scipy.stats import wasserstein_distance

    if not len(angles) == len(adf_ref) == len(adf_pred):
        raise ValueError(f"{len(angles)=}, {len(adf_ref)=}, {len(adf_pred)=} differ")

    angles, adf_ref, adf_pred = map(np.asarray, (angles, adf_ref, adf_pred))
    angle_steps = np.diff(angles)
    if len(angle_steps) == 0:
        raise ValueError("Need >= 2 ADF angle grid points")
    if not np.all(angle_steps > 0) or not np.allclose(angle_steps, angle_steps[0]):
        raise ValueError("ADF angle grid must be strictly increasing and evenly spaced")
    for label, adf in (("reference ADF", adf_ref), ("predicted ADF", adf_pred)):
        _check_density(adf, label)

    # sin(theta) random-neighbor background = the angular ideal gas.
    # wasserstein_distance normalizes weights, so unnormalized densities are fine.
    background = np.sin(np.radians(angles))
    w1_pred = wasserstein_distance(
        angles, angles, u_weights=adf_ref, v_weights=adf_pred
    )
    w1_background = wasserstein_distance(
        angles, angles, u_weights=adf_ref, v_weights=background
    )
    if w1_background <= 0:  # featureless reference -> structural error undefined
        return 100.0
    return float(min(1, w1_pred / w1_background) * 100)


def calc_velocities(trajectory: TrajectoryLike, *, time_step_fs: float) -> np.ndarray:
    """Finite-difference velocities of shape (n_frames, n_atoms, 3) in Angstrom/fs.

    Positions are unwrapped via minimum-image displacements between consecutive
    frames before differentiating, so periodic boundary crossings don't produce
    velocity spikes. Using the same estimator for reference and MLIP trajectories
    keeps vDOS comparisons unbiased even when one stores explicit velocities.
    """
    traj = _as_trajectory(trajectory)
    if traj.n_frames < 2:
        raise ValueError(
            f"Need >= 2 frames to estimate velocities, got {traj.n_frames}"
        )
    if time_step_fs <= 0:
        raise ValueError(f"{time_step_fs=} must be positive")

    positions = traj.positions
    mic_displacements = [
        find_mic(positions[idx + 1] - positions[idx], traj.cell[idx], traj.pbc)[0]
        for idx in range(traj.n_frames - 1)
    ]
    unwrapped = positions[0] + np.cumsum(
        [np.zeros_like(positions[0]), *mic_displacements], axis=0
    )
    return np.gradient(unwrapped, time_step_fs, axis=0)


def equipartition_temperature(
    trajectory: TrajectoryLike, *, time_step_fs: float
) -> float:
    """Mean kinetic temperature (K) implied by finite-difference velocities.

    From equipartition, the mean kinetic energy is ``(3/2) N k_B T``, so
    ``T = mean(sum_i m_i |v_i|^2) / (3 N k_B)``. This is a cheap physical
    sanity check on the saved-frame interval: velocities scale
    as 1/dt, so kinetic energy and hence T scale as 1/dt^2. A reference trajectory whose
    ``dt_fs`` is wrong by a factor f therefore reports a temperature off by f^2 from its
    thermostat target (e.g. dt 5x too large -> T 25x too low), which is how a mislabeled
    ``dt_fs`` silently corrupts the dt-dependent vDOS metric. Pure finite differencing
    slightly underestimates T for trajectories with modes near the Nyquist limit, so
    compare against the target with a generous tolerance, not for exact equality.
    """
    traj = _as_trajectory(trajectory)
    velocities = calc_velocities(traj, time_step_fs=time_step_fs)  # Angstrom/fs
    masses = atomic_masses[traj.atomic_numbers]  # amu
    # 1 amu * (Angstrom/fs)^2 in eV.
    amu_ang2_fs2_in_ev = 1e10 * units.J / units.kg
    kinetic_ev = 0.5 * amu_ang2_fs2_in_ev * np.einsum("a,tad->t", masses, velocities**2)
    return float((2 * kinetic_ev / (3 * traj.n_atoms * units.kB)).mean())


def calc_vdos(
    velocities: np.ndarray, *, time_step_fs: float
) -> tuple[np.ndarray, np.ndarray]:
    """Vibrational density of states from velocities of shape (n_frames, n_atoms, 3).

    Returns frequencies in THz and vDOS intensities: the one-sided power spectrum of
    each Hann-windowed velocity component (the Fourier transform of its
    autocorrelation function by the Wiener-Khinchin theorem), normalized to unit
    total power per component before summing. Per-component normalization weights
    every degree of freedom equally as required by equipartition (a raw sum would
    overweight light atoms by 1/mass) and mirrors the paper's normalized-VACF
    estimator. Components with zero total power (frozen atoms) are skipped.
    """
    velocities = np.asarray(velocities, dtype=float)
    if velocities.ndim != 3 or velocities.shape[-1] != 3:
        raise ValueError(
            f"Expected shape (n_frames, n_atoms, 3), got {velocities.shape}"
        )
    # hanning(n < 4) has at most one non-zero sample, making the windowed
    # spectrum identically zero or a featureless constant
    if (n_frames := len(velocities)) < 4:
        raise ValueError(f"Need >= 4 frames to compute vDOS, got {n_frames}")
    if time_step_fs <= 0:
        raise ValueError(f"{time_step_fs=} must be positive")

    window = np.hanning(n_frames)[:, None, None]
    spectra = np.abs(np.fft.rfft(velocities * window, axis=0)) ** 2
    total_power = spectra.sum(axis=0, keepdims=True)
    spectra = np.divide(
        spectra, total_power, out=np.zeros_like(spectra), where=total_power > 0
    )
    freqs_thz = np.fft.rfftfreq(n_frames, d=time_step_fs * 1e-3)  # fs -> ps = 1/THz
    return freqs_thz, spectra.sum(axis=(1, 2))


def calc_vdos_error(
    freqs_ref: np.ndarray,
    vdos_ref: np.ndarray,
    freqs_pred: np.ndarray,
    vdos_pred: np.ndarray,
    *,
    band_frac: float = 0.999,
) -> float:
    """Wasserstein vDOS error in percent: ``100 * min(1, W1(ref, pred) / sigma_ref)``.
    Matbench W1 variant of the paper's vDOS error (its Eq. 10 uses L1 spectral overlap).

    Wasserstein-1 (earth-mover) distance between the area-normalized reference and MLIP
    vibrational power spectra, normalized by the reference's own spectral std
    ``sigma_ref`` (its frequency spread). Unlike an L1 overlap, W1 grows with how far
    vibrational weight is displaced (systematic phonon softening/hardening) instead of
    saturating (or even inverting) once peaks stop overlapping; on this benchmark it
    roughly doubles the model-discriminating spread (see metrics/md-metrics-design.md).
    The reference high-frequency
    tail beyond ``band_frac`` of its cumulative power is clipped first so W1 is not
    dominated by spectral noise far from the physical bands. 0% = perfect match, 100% =
    displaced by at least the reference's own spectral width (capped). A degenerate
    single-frequency reference (zero spread) returns 0 if the prediction matches it else
    100; a prediction with no power in the reference band also returns 100.
    """
    from scipy.stats import wasserstein_distance

    if not 0 < band_frac <= 1:
        raise ValueError(f"band_frac must be in (0, 1], got {band_frac}")
    f_max = min(np.max(freqs_ref), np.max(freqs_pred))
    grid = np.union1d(freqs_ref, freqs_pred)
    grid = grid[grid <= f_max]
    if len(grid) < 2:
        raise ValueError("Frequency grids have insufficient overlap")

    ref = np.interp(grid, freqs_ref, vdos_ref)
    pred = np.interp(grid, freqs_pred, vdos_pred)
    for label, spectrum in (("reference vDOS", ref), ("predicted vDOS", pred)):
        _check_density(spectrum, label)

    # convert densities to trapezoidal bin masses (density x local cell width) so the
    # score tracks spectral mass, not the (generally non-uniform) merged-grid sampling.
    # Trapezoidal half-cells at the endpoints avoid over-weighting the omega=0 (DC) bin,
    # which carries real weight for diffusive/high-temperature systems.
    bin_widths = np.empty_like(grid)
    bin_widths[1:-1] = (grid[2:] - grid[:-2]) / 2
    bin_widths[[0, -1]] = (grid[1] - grid[0]) / 2, (grid[-1] - grid[-2]) / 2
    ref_mass, pred_mass = ref * bin_widths, pred * bin_widths

    # clip the noisy high-frequency tail beyond band_frac of the reference's power
    # (band_frac == 1 keeps the full grid, even with exact zero-power tail bins)
    if band_frac < 1:
        cut_idx = min(
            int(np.searchsorted(np.cumsum(ref_mass), band_frac * ref_mass.sum())),
            len(grid) - 1,
        )
        mask = grid <= grid[cut_idx]
        grid, ref_mass, pred_mass = grid[mask], ref_mass[mask], pred_mass[mask]
    if ref_mass.sum() <= 0:  # defensive: clipping keeps >= band_frac of reference power
        raise ValueError("reference vDOS has non-positive area after band-clipping")
    if pred_mass.sum() <= 0:  # no MLIP power in the reference band -> max disagreement
        return 100.0

    # wasserstein_distance normalizes weights, so unnormalized masses are fine
    w1 = wasserstein_distance(grid, grid, u_weights=ref_mass, v_weights=pred_mass)
    mean_freq = np.average(grid, weights=ref_mass)
    sigma_ref = float(np.sqrt(np.average((grid - mean_freq) ** 2, weights=ref_mass)))
    if sigma_ref <= 0:  # degenerate single-frequency reference (no spread)
        return 0.0 if w1 <= 0 else 100.0
    return float(min(1, w1 / sigma_ref) * 100)


def get_trajectory_pressures(
    trajectory: TrajectoryLike,
) -> np.ndarray:
    """Per-frame pressures in GPa from the per-frame Voigt stresses (paper Eq. 5),
    vectorized over frames. Raises ValueError if the trajectory has no stress data.
    """
    traj = _as_trajectory(trajectory)
    if traj.stress is None:
        raise ValueError("trajectory has no stress data to compute pressures")
    return -traj.stress[:, :3].sum(axis=1) / 3 / units.GPa


def _time_strides(ref_time_step_fs: float, pred_time_step_fs: float) -> tuple[int, int]:
    """Frame strides (ref, pred) that subsample two trajectories with different
    saved-frame intervals onto their common (least common multiple) time grid, i.e.
    frames i * ref_stride and i * pred_stride share identical timestamps.
    """
    if min(ref_time_step_fs, pred_time_step_fs) <= 0:
        raise ValueError(
            f"Time steps must be positive, got {ref_time_step_fs=}, "
            f"{pred_time_step_fs=}"
        )
    ref_time_step = Fraction(str(ref_time_step_fs))
    pred_time_step = Fraction(str(pred_time_step_fs))
    common_time_step = Fraction(
        math.lcm(ref_time_step.numerator, pred_time_step.numerator),
        math.gcd(ref_time_step.denominator, pred_time_step.denominator),
    )
    return int(common_time_step / ref_time_step), int(common_time_step / pred_time_step)


def calc_pressure_metrics(
    pressures_ref: np.ndarray,
    pressures_pred: np.ndarray,
) -> dict[str, float]:
    """All pressure metrics: the mean-pressure bias (MAE) and Wasserstein-1 distance in
    GPa, plus the histogram-overlap error E_P (paper Eq. 9) in percent.

    The MAE is the absolute difference of the separately-averaged trajectory pressures,
    ``|mean(P_ref) - mean(P_pred)|`` (paper Eq. 6-7). It deliberately does not pair
    frames: reference and MLIP trajectories are independently thermalized, so paired
    instantaneous pressures decorrelate and a frame-paired MAE would measure fluctuation
    noise rather than model bias. W1 and E_P compare the full pressure distributions,
    capturing both the mean offset and fluctuation-width differences.
    """
    from scipy.stats import wasserstein_distance

    pressures_ref, pressures_pred = map(np.asarray, (pressures_ref, pressures_pred))
    if len(pressures_ref) == 0 or len(pressures_pred) == 0:
        raise ValueError("Cannot compute pressure metrics of empty pressure arrays")

    mae = abs(float(pressures_ref.mean()) - float(pressures_pred.mean()))
    return {
        "pressure_mae": mae,
        "pressure_wasserstein": float(
            wasserstein_distance(pressures_ref, pressures_pred)
        ),
        "pressure_error": calc_pressure_histogram_error(pressures_ref, pressures_pred),
    }


def calc_pressure_histogram_error(
    pressures_ref: np.ndarray, pressures_pred: np.ndarray, *, n_bins: int = 80
) -> float:
    """Pressure-distribution error E_P (paper Eq. 9): the non-overlap of the
    area-normalized reference and MLIP pressure histograms over shared bin edges, in
    percent. 0% = identical distributions, 100% = disjoint.

    Eq. 9's denominator is identically 2 (each density integrates to 1), so this is
    ``50 * sum |density_ref - density_pred| * bin_width``.
    """
    pressures_ref, pressures_pred = map(np.asarray, (pressures_ref, pressures_pred))
    if n_bins < 2:
        raise ValueError(f"Need >= 2 pressure histogram bins, got {n_bins=}")
    if len(pressures_ref) == 0 or len(pressures_pred) == 0:
        raise ValueError("Cannot compute pressure histogram error of empty arrays")
    all_pressures = np.concatenate([pressures_ref, pressures_pred])
    if not np.all(np.isfinite(all_pressures)):
        raise ValueError("Pressure histogram error requires finite pressure arrays")
    lo, hi = float(all_pressures.min()), float(all_pressures.max())
    if hi == lo:  # both distributions collapse to one identical value -> full overlap
        return 0.0
    edges = np.linspace(lo, hi, n_bins + 1)
    density_ref, _ = np.histogram(pressures_ref, bins=edges, density=True)
    density_pred, _ = np.histogram(pressures_pred, bins=edges, density=True)
    return float(50 * np.sum(np.abs(density_ref - density_pred) * np.diff(edges)))


def calc_energy_force_rmse(
    trajectory: TrajectoryLike, calculator: "Calculator"
) -> dict[str, float]:
    """Private-label energy-fluctuation (eV/atom) and force (eV/Å) RMSE of calculator
    single-points on labeled reference frames.

    The energy RMSE uses per-trajectory mean-subtracted energies, i.e. the energy
    *fluctuations* sampled along the trajectory, not absolute energies (invariant to
    the reference's arbitrary energy zero). Forces are reference-zero-invariant, so
    their RMSE stays absolute.
    """
    traj = _as_trajectory(trajectory)
    if traj.n_frames == 0:
        raise ValueError("Cannot evaluate a calculator on an empty trajectory")
    if traj.energy is None or traj.forces is None:
        raise ValueError("trajectory lacks private reference energy/forces for RMSE")

    energy_pred = np.empty(traj.n_frames)
    force_se = np.empty(traj.n_frames)
    for frame_idx in range(traj.n_frames):
        # fresh Atoms (no stored results) so the calculator computes the prediction
        atoms = Atoms(
            numbers=traj.atomic_numbers,
            positions=traj.positions[frame_idx],
            cell=traj.cell[frame_idx],
            pbc=traj.pbc,
        )
        atoms.calc = calculator
        energy_pred[frame_idx] = atoms.get_potential_energy()
        force_se[frame_idx] = ((atoms.get_forces() - traj.forces[frame_idx]) ** 2).sum()

    pred_dev = (energy_pred - energy_pred.mean()) / traj.n_atoms
    ref_dev = (traj.energy - traj.energy.mean()) / traj.n_atoms
    n_force_components = traj.n_frames * traj.n_atoms * 3
    return {
        "energy_rmse": float(np.sqrt(np.mean((pred_dev - ref_dev) ** 2))),
        "force_rmse": float(np.sqrt(force_se.sum() / n_force_components)),
    }


def matched_frame_counts(
    *,
    n_ref_frames: int,
    n_pred_frames: int,
    ref_time_step_fs: float,
    pred_time_step_fs: float,
) -> tuple[int, int]:
    """Frame counts (n_ref_use, n_pred_use) covering the identical closed time
    interval [0, t_match] in both trajectories, where t_match is the largest
    multiple of the common (LCM) time grid that fits both trajectory durations.
    A trajectory of n frames saved every dt spans (n - 1) * dt.
    """
    ref_stride, pred_stride = _time_strides(ref_time_step_fs, pred_time_step_fs)
    if min(n_ref_frames, n_pred_frames) < 1:
        return 0, 0
    # number of common-grid steps fitting both trajectory durations
    n_common = min((n_ref_frames - 1) // ref_stride, (n_pred_frames - 1) // pred_stride)
    # + 1 to include the frame at t=0
    return n_common * ref_stride + 1, n_common * pred_stride + 1


def evaluate_md_system(
    ref_trajectory: TrajectoryLike,
    pred_trajectory: TrajectoryLike,
    *,
    ref_time_step_fs: float,
    pred_time_step_fs: float,
    n_rdf_bins: int = 500,
    progress_label: str | None = None,
) -> dict[str, float]:
    """All per-system MD metrics for one reference/predicted trajectory pair, keyed
    by METRIC_UNITS column names.

    Both trajectories are truncated to the same total simulation time (given their
    respective time steps in fs between saved frames) before computing observables,
    matching the paper's same-simulation-length protocol. Pressure metrics are NaN
    when either trajectory lacks stress data.

    Args:
        ref_trajectory: Ab-initio reference trajectory.
        pred_trajectory: Model-predicted trajectory.
        ref_time_step_fs: Time between saved reference frames in fs.
        pred_time_step_fs: Time between saved predicted frames in fs.
        n_rdf_bins: Number of RDF histogram bins.
        progress_label: Optional label for coarse per-metric progress logging.
    """

    def log_progress(stage: str) -> None:
        """Print a coarse metric-stage marker for long Slurm jobs."""
        if progress_label:
            print(f"{progress_label}: {stage}", flush=True)

    ref = _as_trajectory(ref_trajectory)
    pred = _as_trajectory(pred_trajectory)
    n_ref_use, n_pred_use = matched_frame_counts(
        n_ref_frames=ref.n_frames,
        n_pred_frames=pred.n_frames,
        ref_time_step_fs=ref_time_step_fs,
        pred_time_step_fs=pred_time_step_fs,
    )
    if n_ref_use < 4 or n_pred_use < 4:  # vDOS needs >= 4 frames
        raise ValueError(
            f"Trajectories too short after time matching: {n_ref_use=}, {n_pred_use=}"
        )
    ref_matched = ref[:n_ref_use]
    pred_matched = pred[:n_pred_use]

    log_progress("RDF")
    radii, g_r_ref = calc_rdf(ref_matched, n_bins=n_rdf_bins)
    r_max = float(radii[-1] + (radii[1] - radii[0]) / 2)
    _, g_r_pred = calc_rdf(pred_matched, n_bins=n_rdf_bins, r_max=r_max)
    metrics = {"rdf_error": calc_rdf_error(radii, g_r_ref, g_r_pred)}
    log_progress("ADF")
    angles, adf_ref = calc_adf(ref_matched)
    _, adf_pred = calc_adf(pred_matched)
    metrics["adf_error"] = calc_adf_error(angles, adf_ref, adf_pred)

    log_progress("vDOS")
    ref_velocities = calc_velocities(ref_matched, time_step_fs=ref_time_step_fs)
    pred_velocities = calc_velocities(pred_matched, time_step_fs=pred_time_step_fs)
    freqs_ref, vdos_ref = calc_vdos(ref_velocities, time_step_fs=ref_time_step_fs)
    freqs_pred, vdos_pred = calc_vdos(pred_velocities, time_step_fs=pred_time_step_fs)
    metrics["vdos_error"] = calc_vdos_error(freqs_ref, vdos_ref, freqs_pred, vdos_pred)

    log_progress("pressure")
    if ref_matched.stress is not None and pred_matched.stress is not None:
        metrics |= calc_pressure_metrics(
            get_trajectory_pressures(ref_matched),
            get_trajectory_pressures(pred_matched),
        )
    else:  # one or both trajectories lack stress -> undefined pressure metrics
        metrics |= dict.fromkeys(
            ("pressure_mae", "pressure_wasserstein", "pressure_error"), float("nan")
        )

    log_progress("metrics done")
    return metrics


def combine_per_system_metrics(frames: Sequence[pd.DataFrame]) -> pd.DataFrame:
    """Concatenate already-read per-system MD metric rows (one row each, as written by
    parallel single-system runs) into a single DataFrame indexed by system.

    Takes pre-read DataFrames rather than paths so callers that already read each CSV
    (e.g. to filter out multi-system files) don't read it twice. Deduplicates on the
    system column keeping the last occurrence, so re-running a system (e.g. a timed-out
    rollout finished in a later job) overrides the earlier row, not double-counting it.
    """
    if not frames:
        raise ValueError("No per-system metric frames given")
    df_all = pd.concat(frames, ignore_index=True)
    if "system" not in df_all:
        raise ValueError(f"per-system metrics lack a 'system' column, got {[*df_all]}")
    return df_all.drop_duplicates(subset="system", keep="last").set_index("system")


def calc_md_metrics(df_md: pd.DataFrame) -> dict[str, float | str]:
    """Aggregate per-system MD metric rows (one row per system with a subset of
    PER_SYSTEM_METRIC_COLS columns) into model-level metrics.

    Means are taken over systems, skipping NaNs (e.g. systems without stress data).
    The combined MD score (CMDS) is deliberately absent: it is site-computed (see
    module docstring). Private-label energy/force RMSE rows are in eV but reported
    here in meV for readability.

    Run provenance columns (RUN_PROVENANCE_COLS) aggregate differently:
    ``run_time_sec`` is summed over systems (parallel per-system jobs ~ one serial
    sweep, mirroring the diatomics shard merge), the memory peaks take the max over
    systems (the largest system sets the hardware requirement) and ``hardware``
    passes through. All are all-or-nothing like the private-label diagnostics: a
    partial sum/max or mixed-GPU label would misrepresent the model's cost, so they
    are only published when every system row agrees.
    """
    metric_cols = [col for col in PER_SYSTEM_METRIC_COLS if col in df_md]
    if not metric_cols:
        raise ValueError(
            f"No recognized MD metric columns in {list(df_md.columns)=}, "
            f"expected a subset of {PER_SYSTEM_METRIC_COLS}"
        )
    metrics: dict[str, float | str] = {}
    # hardware is uniform when every row has the same non-null label; an incomplete or
    # mixed column also drops the GPU-dependent costs (summed wall time, peak VRAM): a
    # total over an H200 and an A100 (or rows of unknown provenance) would misrepresent
    # the model's cost with no hardware label left to qualify it. Host RSS is kept: it
    # reflects model/system size, not the GPU.
    hardware_uniform = True
    if "hardware" in df_md:
        hardware_vals = df_md["hardware"].unique()
        hardware_uniform = len(hardware_vals) == 1 and pd.notna(hardware_vals[0])
        if hardware_uniform:
            metrics["hardware"] = str(hardware_vals[0])
    for col, agg, gpu_dependent in (
        ("run_time_sec", "sum", True),
        ("max_rss_gb", "max", False),
        ("max_gpu_mem_gb", "max", True),
    ):
        if gpu_dependent and not hardware_uniform:
            continue
        if col in df_md and df_md[col].notna().all():
            metrics[col] = float(df_md[col].agg(agg))
    metrics |= {
        col: float(df_md[col].mean())
        * (1000 if col in ("energy_rmse", "force_rmse") else 1)
        for col in metric_cols
    }
    # private-label diagnostics must cover every system or not be published at all: a
    # NaN-skipping subset mean (mixed runs with/without a private ref) would silently
    # misrepresent the model while n_systems suggests full coverage. (Pressure NaNs
    # differ: stress-less systems are NaN for every model alike, means stay comparable.)
    for col in ("energy_rmse", "force_rmse"):
        if col in metrics and df_md[col].isna().any():
            del metrics[col]
    metrics["n_systems"] = len(df_md)
    return metrics


def write_metrics_to_yaml(
    model: "Model",
    metrics: dict[str, float | str],
    pred_file_path: str | None = None,
    pred_file_url: str | None = None,
) -> dict[str, Any]:
    """Write model-level MD metrics from calc_md_metrics to the metrics.md section
    of a model YAML file, with unit comments and optional pred_file references.
    """
    yaml_metrics = CommentedMap()
    if pred_file_path is not None:
        yaml_metrics["pred_file"] = repo_relative_path(pred_file_path)
    if pred_file_url is not None:
        yaml_metrics["pred_file_url"] = pred_file_url
    for key, value in metrics.items():
        if key == "n_systems":
            yaml_metrics[key] = int(value)
        elif isinstance(value, str):  # e.g. hardware
            yaml_metrics[key] = value
        else:
            yaml_metrics[key] = float(round(value, 4))
        if unit := METRIC_UNITS.get(key):
            yaml_metrics.yaml_add_eol_comment(unit, key, column=1)

    update_yaml_file(model.yaml_path, "metrics.md", yaml_metrics)
    return yaml_metrics
