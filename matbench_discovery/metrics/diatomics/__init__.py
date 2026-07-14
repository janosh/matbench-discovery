"""Metrics for analyzing diatomic potential energy curves.

Big props to Tamas Stenczel and Yuan Chiang who spearheaded this type of PES
smoothness analysis in https://github.com/stenczelt/MACE-MP-work for the MACE-MP
paper https://arxiv.org/abs/2401.00096 (see fig. 56) and MLIP Arena
https://huggingface.co/spaces/atomind/mlip-arena, respectively.
"""

import gzip
import json
from dataclasses import dataclass, field
from typing import Any, Self

import numpy as np
from ase.data import atomic_numbers, covalent_radii, vdw_alvarez
from numpy.typing import ArrayLike

from matbench_discovery import repo_relative_path
from matbench_discovery.data import (
    FileRef,
    file_ref_name,
    file_ref_url,
    make_file_ref,
    update_yaml_file,
)
from matbench_discovery.enums import MbdKey, Model
from matbench_discovery.metrics.diatomics.energy import (
    calc_energy_diff_flips,
    calc_energy_jump,
    calc_pbe_bond_length_error,
    calc_pbe_energy_mae,
    calc_pbe_vib_freq_error,
    calc_pbe_wall_dist_mae,
    calc_pbe_well_depth_error,
    calc_tortuosity,
)
from matbench_discovery.metrics.diatomics.force import (
    calc_force_flips,
    calc_force_jump,
    calc_force_mae,
    calc_force_total_variation,
)

DiatomicsYamlValue = str | float | dict[str, str] | FileRef | None
# Elements absent from the Materials Project (MP covers 89 elements: H-Pu minus these
# five), and hence from MPtrj/OMat24-derived training data. Models trained on MP data
# cannot predict them, so diatomic metrics skip them for every model (and the runner
# tolerates their missing curves without per-model exclusion bookkeeping).
NON_MP_ELEMENTS = frozenset({"Po", "At", "Rn", "Fr", "Ra"})
DIATOMIC_WALL_R_MIN_FACTOR = 0.8
DIATOMIC_METRIC_KEYS = frozenset(
    str(key)
    for key in (
        MbdKey.tortuosity,
        MbdKey.force_flips,
        MbdKey.energy_jump,
        MbdKey.energy_diff_flips,
        MbdKey.force_total_variation,
        MbdKey.force_jump,
        MbdKey.pbe_wall_dist_mae,
        MbdKey.pbe_energy_mae,
        MbdKey.pbe_bond_length_error,
        MbdKey.pbe_well_depth_error,
        MbdKey.pbe_force_mae,
        MbdKey.pbe_vib_freq_error,
    )
)


def _homo_key(formula: str) -> str:
    """Return element key for homonuclear pair labels, else formula unchanged."""
    elem1, sep, elem2 = formula.partition("-")
    return elem1 if sep and elem1 == elem2 else formula


class DiatomicCurve:
    """Energies and forces for a single diatomic molecule at multiple distances."""

    distances: np.ndarray  # shape (n_distances,)
    energies: np.ndarray  # shape (n_distances,)
    forces: np.ndarray  # shape (n_distances, n_atoms, 3)

    def __init__(
        self, distances: ArrayLike, energies: ArrayLike, forces: ArrayLike
    ) -> None:
        """Convert inputs to numpy arrays and reshape forces if needed."""
        self.distances = np.asarray(distances)
        self.energies = np.asarray(energies)
        self.forces = np.asarray(forces)

        n_dist, n_e = len(self.distances), len(self.energies)
        if n_e != n_dist:
            raise ValueError(f"distance and energy counts differ: {n_dist} != {n_e}")

        # Handle forces stored as (1, n_distances*n_atoms, 3)
        # instead of (n_distances, n_atoms, 3)
        if self.forces.shape == (1, 2 * n_dist, 3):
            self.forces = self.forces.reshape(n_dist, 2, 3)
        if (n_f := len(self.forces)) != n_dist:
            raise ValueError(f"distance and force counts differ: {n_dist} != {n_f}")


@dataclass
class DiatomicCurves:
    """Container for diatomic potential energy curves and forces of multiple
    element pairs.

    Attributes:
        distances (np.ndarray): Interatomic distances in Å.
        homo_nuclear (dict[str, DiatomicCurve]): Map of element pairs
            (e.g. "H-H") to their DiatomicCurve (energies and forces).
        hetero_nuclear (dict[str, DiatomicCurve] | None): Optional map of element pairs
            (e.g. "H-He") to their DiatomicCurve.
    """

    distances: np.ndarray  # shape (n_distances,)
    homo_nuclear: dict[str, DiatomicCurve]
    hetero_nuclear: dict[str, DiatomicCurve] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        """Create DiatomicCurves from a dictionary loaded from JSON."""
        distances = np.asarray(data["distances"])
        grid_pos_by_dist = {float(dist): idx for idx, dist in enumerate(distances)}

        def make_curves(section: str) -> dict[str, DiatomicCurve]:
            """Convert one JSON curve section to DiatomicCurve objects."""
            curves = data.get(section, data.get(section.replace("-", "_"), {}))
            key_fn = _homo_key if section.startswith("homo") else str

            def curve_distances(formula: str, curve: dict[str, Any]) -> np.ndarray:
                """Return curve distances, requiring consistency with top-level grid.

                Curves may omit grid points (run_diatomics trims non-finite samples
                outside the scored window) but must not introduce off-grid,
                duplicated, or reordered distances.
                """
                curve_dists = np.asarray(curve.get("distances", distances))
                # off-grid points map to -1; valid subsets have strictly
                # increasing grid positions
                grid_positions = np.array(
                    [grid_pos_by_dist.get(float(dist), -1) for dist in curve_dists]
                )
                if (grid_positions < 0).any() or (np.diff(grid_positions) <= 0).any():
                    raise ValueError(
                        f"{formula} curve distances must be an ordered subset of "
                        "top-level distances"
                    )
                return curve_dists

            return {
                key_fn(formula): DiatomicCurve(
                    distances=curve_distances(formula, curve),
                    energies=curve["energies"],
                    forces=curve.get("forces", []),
                )
                for formula, curve in curves.items()
                if len(curve["energies"]) > 0
            }

        return cls(
            distances=distances,
            homo_nuclear=make_curves("homo-nuclear"),
            hetero_nuclear=make_curves("hetero-nuclear"),
        )


def load_dft_reference_curves(
    functional: str = "PBE", ref_path: str | None = None
) -> DiatomicCurves:
    """Load bundled DFT reference curves as DiatomicCurves keyed by formula."""
    from matbench_discovery import ROOT

    dft_ref_path = ref_path or f"{ROOT}/site/src/lib/diatomics-dft.json.gz"
    open_fn = gzip.open if dft_ref_path.endswith(".gz") else open
    with open_fn(dft_ref_path, mode="rt", encoding="utf-8") as file:
        references = json.load(file)[functional]
    return DiatomicCurves(
        distances=np.array([]),
        homo_nuclear={
            _homo_key(formula): DiatomicCurve(
                distances=curve["distances"],
                energies=curve["energies"],
                forces=curve.get("forces", []),
            )
            for formula, curve in references.items()
        },
    )


def find_low_quality_dft_refs(
    ref_curves: DiatomicCurves,
    *,
    min_energy_jump: float = 1.5,
    min_energy_flips: int = 3,
) -> set[str]:
    """Elements whose DFT reference curve is too jumpy to score models against.

    A reference curve is flagged when its own smoothness metrics (evaluated over the
    same element-specific window used to score models) mark it as oscillatory: total
    energy jump at sign-flip points >= min_energy_jump AND number of energy-difference
    sign flips >= min_energy_flips. Requiring both avoids flagging curves with a few
    large but possibly physical features (e.g. Cr2's shelf) or many tiny numerical
    wiggles. On the current PBE reference, this flags 8 lanthanides (Pr, Pm, Sm,
    Tb, Dy, Ho, Er, Tm) whose f-electron SCF convergence issues produce eV-scale
    discontinuities, next to which any model error signal drowns.

    Args:
        ref_curves (DiatomicCurves): DFT reference curves keyed by element symbol.
        min_energy_jump (float): Flag threshold for calc_energy_jump in eV.
        min_energy_flips (int): Flag threshold for calc_energy_diff_flips.

    Returns:
        set[str]: Element symbols with low-quality reference curves.
    """
    low_quality: set[str] = set()
    for elem_symbol, curve in ref_curves.homo_nuclear.items():
        seps, energies = curve.distances, curve.energies  # np arrays per DiatomicCurve
        if seps.size == 0:
            continue
        r_min, r_max = eval_window(elem_symbol, float(np.max(seps)))
        mask = (seps >= r_min) & (seps <= r_max)
        if mask.sum() < 5:
            continue  # too few in-window points to assess smoothness
        if not np.isfinite(energies[mask]).all():
            # non-finite reference energies can't be scored against (and would crash
            # the smoothness metrics below)
            low_quality.add(elem_symbol)
            continue
        if (
            calc_energy_jump(seps[mask], energies[mask]) >= min_energy_jump
            and calc_energy_diff_flips(seps[mask], energies[mask]) >= min_energy_flips
        ):
            low_quality.add(elem_symbol)
    return low_quality


def eval_window(
    elem_symbol: str, seps_max: float, *, r_min_factor: float = 0.9
) -> tuple[float, float]:
    """Return the element-specific metric-evaluation window.

    The default lower bound matches MLIP Arena at ``0.9 * covalent_radius`` so the
    steep repulsive wall does not dominate every magnitude metric. Reference-relative
    wall scoring passes ``r_min_factor=0.8`` to include the full DFT-supported range.

    Args:
        elem_symbol (str): Homonuclear pair label, e.g. "H-H".
        seps_max (float): Largest available separation, used to cap r_max.
        r_min_factor (float): Lower-bound multiplier on the covalent radius.

    Returns:
        tuple[float, float]: (r_min, r_max) in Å.
    """
    atomic_num = atomic_numbers[elem_symbol.split("-", maxsplit=1)[0]]
    r_cov = covalent_radii[atomic_num] if atomic_num < len(covalent_radii) else np.nan
    r_min = r_min_factor * r_cov if np.isfinite(r_cov) else 0.0
    vdw_radii = vdw_alvarez.vdw_radii
    r_vdw = vdw_radii[atomic_num] if atomic_num < len(vdw_radii) else np.nan
    r_max = min(3.1 * r_vdw, seps_max) if np.isfinite(r_vdw) else seps_max
    return r_min, r_max


def calc_diatomic_metrics(
    ref_curves: DiatomicCurves | None,
    pred_curves: DiatomicCurves,
    metrics: dict[str, dict[str, Any]] | None = None,
    *,
    interpolate: bool | int = False,
) -> dict[str, dict[str, float]]:
    """Calculate diatomic curve metrics comparing predicted curves to reference curves.

    Args:
        ref_curves (DiatomicCurves | None): Reference energy curves for each element.
            If None, only metrics that don't require reference data will be calculated.
        pred_curves (DiatomicCurves): Predicted energy curves for each element.
        metrics (dict[str, dict[str, Any]] | None): Map of metric names to
            dictionaries of keyword arguments for each metric function. If None, uses
            all metrics with default parameters. To use a subset of metrics, provide
            a dictionary with those metric names as keys and their keyword arguments
            as values. Empty dictionaries will use default parameters.
        interpolate (bool | int): If False (default), uses the provided points directly.
            If True, uses default number of points for interpolation.
            If an integer, uses that many points for interpolation.
            This value is passed to each metric function that supports it.

    Returns:
        dict[str, dict[str, float]]: Map of element symbols to metric dicts with keys
            being the metric names and values being the metric values. Elements outside
            the Materials Project set (NON_MP_ELEMENTS) are skipped entirely since
            MP-trained models cannot support them. Elements whose DFT reference curve
            is too jumpy to score against (see find_low_quality_dft_refs) still get
            self-consistency metrics but no reference-relative pbe_* metrics.
    """
    requested_metric_keys = (
        {str(key) for key in metrics} if metrics is not None else DIATOMIC_METRIC_KEYS
    )
    if unknown_metrics := requested_metric_keys - DIATOMIC_METRIC_KEYS:
        raise ValueError(f"{unknown_metrics=}. Valid metrics=")
    requested_metrics = {MbdKey(key) for key in requested_metric_keys}

    results: dict[str, dict[str, float]] = {}
    metric_kwargs = {
        MbdKey(key): kwargs.copy() for key, kwargs in (metrics or {}).items()
    }
    for key in requested_metrics & {
        MbdKey.pbe_energy_mae,
        MbdKey.pbe_force_mae,
    }:
        metric_kwargs.setdefault(key, {}).setdefault("interpolate", interpolate)

    low_quality_refs = find_low_quality_dft_refs(ref_curves) if ref_curves else set()
    for elem_symbol, pred_data in pred_curves.homo_nuclear.items():
        if _homo_key(elem_symbol) in NON_MP_ELEMENTS:
            continue  # score all models on the same MP-supported element set
        # Restrict general metrics to MLIP Arena's element-specific window so the steep
        # repulsive wall does not dominate every magnitude metric. The dedicated wall
        # metric below uses the full DFT range down to 0.8 * covalent radius.
        pred_dists = pred_data.distances
        seps_max = float(pred_dists.max())
        r_min, r_max = eval_window(elem_symbol, seps_max)
        pred_mask = (pred_dists >= r_min) & (pred_dists <= r_max)
        if pred_mask.sum() < 5:  # too few points in window for stable metrics
            print(f"Skipping {elem_symbol} diatomic metrics: <5 points in eval window")
            continue
        seps_pred = pred_dists[pred_mask]
        pred_energies_raw = np.asarray(pred_data.energies)
        pred_energies = pred_energies_raw[pred_mask]
        pred_forces_raw = np.asarray(pred_data.forces)
        if not pred_forces_raw.size:
            raise ValueError(f"{elem_symbol} diatomic curve is missing forces")
        if len(pred_forces_raw) != len(pred_dists):
            raise ValueError(
                f"{elem_symbol} diatomic force and distance counts differ: "
                f"{len(pred_forces_raw)} != {len(pred_dists)}"
            )
        pred_forces = pred_forces_raw[pred_mask]

        # skip curves still non-finite within the window (model instability at scored
        # geometries); they'd fail metric validation and abort the whole model otherwise
        if not (np.isfinite(pred_energies).all() and np.isfinite(pred_forces).all()):
            print(f"Skipping {elem_symbol} diatomic metrics: non-finite curve values")
            continue

        energy_args = (seps_pred, pred_energies)
        force_args = (seps_pred, pred_forces)
        # (metric_key, function, positional args) for metrics needing only the
        # predicted curve
        metric_calls: list[tuple[MbdKey, Any, tuple[Any, ...]]] = [
            (MbdKey.tortuosity, calc_tortuosity, energy_args),
            (MbdKey.energy_diff_flips, calc_energy_diff_flips, energy_args),
            (MbdKey.energy_jump, calc_energy_jump, energy_args),
            (MbdKey.force_flips, calc_force_flips, force_args),
            (MbdKey.force_total_variation, calc_force_total_variation, force_args),
            (MbdKey.force_jump, calc_force_jump, force_args),
        ]

        # prepend reference-requiring metrics when a matching reference curve exists
        # and passes the smoothness quality gate (jumpy references would drown any
        # model error signal, but self-consistency metrics above stay meaningful)
        if (
            ref_curves
            and _homo_key(elem_symbol) not in low_quality_refs
            and (ref_data := ref_curves.homo_nuclear.get(elem_symbol))
        ):
            ref_dists = ref_data.distances
            ref_mask = (ref_dists >= r_min) & (ref_dists <= r_max)
            seps_ref = ref_dists[ref_mask]
            ref_energies_raw = np.asarray(ref_data.energies)
            ref_energies = ref_energies_raw[ref_mask]
            if len(seps_ref) >= 2:
                pair_args = (seps_ref, ref_energies, seps_pred, pred_energies)
                metric_calls[:0] = [
                    (MbdKey.pbe_energy_mae, calc_pbe_energy_mae, pair_args),
                    (
                        MbdKey.pbe_bond_length_error,
                        calc_pbe_bond_length_error,
                        pair_args,
                    ),
                    (MbdKey.pbe_well_depth_error, calc_pbe_well_depth_error, pair_args),
                    (
                        MbdKey.pbe_vib_freq_error,
                        calc_pbe_vib_freq_error,
                        (elem_symbol, *pair_args),
                    ),
                ]
            wall_r_min = eval_window(
                elem_symbol, seps_max, r_min_factor=DIATOMIC_WALL_R_MIN_FACTOR
            )[0]
            # The generated DFT endpoint can differ from exactly 0.8*r_cov by one
            # floating-point ulp, so include numerically equal boundary points.
            wall_r_min -= 1e-12
            pred_wall_mask = (pred_dists >= wall_r_min) & (pred_dists <= r_max)
            ref_wall_mask = (ref_dists >= wall_r_min) & (ref_dists <= r_max)
            if pred_wall_mask.sum() >= 2 and ref_wall_mask.sum() >= 2:
                wall_args = (
                    ref_dists[ref_wall_mask],
                    ref_energies_raw[ref_wall_mask],
                    pred_dists[pred_wall_mask],
                    pred_energies_raw[pred_wall_mask],
                )
                metric_calls.insert(
                    0, (MbdKey.pbe_wall_dist_mae, calc_pbe_wall_dist_mae, wall_args)
                )
            ref_forces = np.asarray(ref_data.forces)
            if (
                ref_forces.size
                and len(ref_forces) == len(ref_dists)
                and len(seps_ref) >= 2
            ):
                ref_forces = ref_forces[ref_mask]
                force_interpolate = metric_kwargs.get(MbdKey.pbe_force_mae, {}).get(
                    "interpolate", False
                )
                same_grid = np.array_equal(seps_ref, seps_pred)
                has_overlap = max(seps_ref.min(), seps_pred.min()) <= min(
                    seps_ref.max(), seps_pred.max()
                )
                if same_grid or (force_interpolate and has_overlap):
                    force_pair_args = (seps_ref, ref_forces, seps_pred, pred_forces)
                    metric_calls.append(
                        (MbdKey.pbe_force_mae, calc_force_mae, force_pair_args)
                    )

        results[elem_symbol] = {
            key: calc_fn(*args, **metric_kwargs.get(key, {}))
            for key, calc_fn, args in metric_calls
            if key in requested_metrics
        }

    return results


def write_metrics_to_yaml(
    model: Model,
    metrics: dict[str, dict[str, float]],
    pred_file_path: str | None = None,
    run_metadata: dict[str, str | float | dict[str, str]] | None = None,
) -> dict[str, DiatomicsYamlValue]:
    """Write diatomic metrics to model YAML file.

    Args:
        model (Model): Model to write metrics for.
        metrics (dict[str, dict[str, float]]): Map of element symbols to dicts of
            metric values.
        pred_file_path (str | None): If given, record this path as
            metrics.diatomics.pred_file. Absolute paths must be inside the repo and are
            converted to repo-relative paths. Otherwise an existing pred_file is
            preserved.
        run_metadata (dict[str, str | float | dict[str, str]] | None): Extra run
            fields (e.g. hardware, run_time_sec, pred_file_url). ``pred_file_url``
            overrides or supplies the prediction file URL. A recompute without
            run_metadata preserves existing values.

    Returns:
        dict[str, DiatomicsYamlValue]: The metrics.diatomics block written (file refs
            and run metadata first, then metric means across all elements).
    """
    # mean of each metric over the elements that have a finite value: skips elements
    # whose windowed curve is degenerate (e.g. tortuosity is NaN for a flat curve), and
    # drops a metric entirely if no element has a finite value rather than writing an
    # invalid `.nan`. Union the keys since elements can differ (only some have a ref).
    mean_metrics: dict[str, DiatomicsYamlValue] = {}
    for metric in dict.fromkeys(
        key for elem_metrics in metrics.values() for key in elem_metrics
    ):
        finite = [
            val
            for elem_metrics in metrics.values()
            if (val := elem_metrics.get(metric)) is not None and np.isfinite(val)
        ]
        if finite:
            mean_metrics[str(metric)] = float(f"{np.mean(finite):.4}")

    # carry over only recognized run metadata (it describes the source run, not the
    # computed metrics, so it stays valid on recalculation)
    existing = model.metrics.get("diatomics", {})
    existing = existing if isinstance(existing, dict) else {}
    run_metadata = run_metadata or {}
    block: dict[str, DiatomicsYamlValue] = {}
    existing_pred_file = existing.get("pred_file")
    existing_name = file_ref_name(existing_pred_file)
    pred_name = existing_name
    if pred_file_path is not None:
        pred_name = repo_relative_path(pred_file_path)
    pred_file_url = run_metadata.get("pred_file_url")
    if pred_file_url is None and pred_name == existing_name:
        pred_file_url = file_ref_url(existing_pred_file)
    if pred_name is not None:
        url = pred_file_url if isinstance(pred_file_url, str) else None
        size = md5 = None
        if pred_name == existing_name and isinstance(existing_pred_file, dict):
            size_val, md5_val = (
                existing_pred_file.get("size"),
                existing_pred_file.get("md5"),
            )
            if isinstance(size_val, int) and isinstance(md5_val, str):
                size, md5 = size_val, md5_val
        block["pred_file"] = make_file_ref(pred_name, url=url, size=size, md5=md5)

    run_metadata_keys = (
        *("hardware", "run_time_sec", "max_rss_gb", "max_gpu_mem_gb"),
        "excluded_formula_reasons",
    )
    for key in run_metadata_keys:
        if (val := run_metadata.get(key)) is None:
            val = existing.get(key)
        if val is not None:
            block[key] = val
    block |= mean_metrics

    # preserve_existing=False so a recompute fully replaces the block, dropping
    # deprecated metrics left in the YAML. This still runs when no finite metrics were
    # produced, preventing stale metric values from surviving.
    update_yaml_file(
        model.yaml_path, "metrics.diatomics", block, preserve_existing=False
    )
    print(f"Wrote {model.label} diatomic metrics to {model.yaml_path}")
    return block
