"""Build and postprocess DFT diatomic reference curves from spin candidates."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Literal

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

SpinCandidate = int | str
CurvePoint = dict[str, Any]


@dataclass(frozen=True)
class CurvePostprocessEdit:
    """Metadata for one DFT reference curve postprocess edit."""

    kind: Literal[
        "isolated_spin_branch_drop",
        "isolated_energy_bump",
        "collapsed_scf_point",
    ]
    distance: float
    original_spin: str
    replacement_spin: str | None
    original_energy: float
    replacement_energy: float | None
    expected_energy: float


def point_energy(point: Mapping[str, Any]) -> float:
    """Return the final, sigma->0 energy from a curve point."""
    energies = point.get("energies")
    if not energies:
        msg = f"curve point lacks energies: {point}"
        raise ValueError(msg)
    energy = float(energies[-1])
    if not math.isfinite(energy):
        msg = f"curve point has non-finite final energy: {energy}"
        raise ValueError(msg)
    return energy


def point_distance(point: Mapping[str, Any]) -> float:
    """Return the finite distance from a curve point."""
    distance = float(point["distance"])
    if not math.isfinite(distance):
        msg = f"curve point has non-finite distance: {distance}"
        raise ValueError(msg)
    return distance


def distance_key(distance: float, *, ndigits: int = 4) -> float:
    """Round a distance to the canonical key used to align spin candidates."""
    return round(distance, ndigits)


def finite_energy_points(points: Sequence[Mapping[str, Any]]) -> list[CurvePoint]:
    """Return points with finite final energies and distances."""
    finite_points: list[CurvePoint] = []
    for point in points:
        try:
            point_energy(point)
            point_distance(point)
        except (KeyError, TypeError, ValueError):
            continue
        finite_points.append(dict(point))
    return finite_points


def drop_collapsed_scf_points(
    points: Sequence[Mapping[str, Any]],
    *,
    min_collapse_ev: float = 20.0,
    spin_candidate: str = "",
) -> tuple[list[CurvePoint], list[CurvePostprocessEdit]]:
    """Drop variationally collapsed SCF points from one candidate curve.

    Some fixed-spin runs (mostly f-electron heavy lanthanides) collapse to spurious
    solutions tens to thousands of eV below the rest of their own branch (e.g. one
    Gd/r2SCAN AFM point at -261 eV between -74 eV neighbors). No smooth PEC has a
    single-point dip of min_collapse_ev; wells are eV-scale. Such points would always
    win the per-distance min-merge, so they must be removed per candidate before
    merging. Iterates so consecutive collapsed points peel off one by one.
    """
    kept = sorted(finite_energy_points(points), key=point_distance)
    edits: list[CurvePostprocessEdit] = []
    while len(kept) > 1:
        for idx, point in enumerate(kept):
            expected_energy = min(
                point_energy(kept[jdx])
                for jdx in (idx - 1, idx + 1)
                if 0 <= jdx < len(kept)
            )
            energy = point_energy(point)
            if energy < expected_energy - min_collapse_ev:
                edits.append(
                    CurvePostprocessEdit(
                        kind="collapsed_scf_point",
                        distance=point_distance(point),
                        original_spin=spin_candidate,
                        replacement_spin=None,
                        original_energy=energy,
                        replacement_energy=None,
                        expected_energy=expected_energy,
                    )
                )
                del kept[idx]
                break
        else:
            break
    return kept, edits


def merge_min_energy_curve(
    candidate_points: Mapping[SpinCandidate, Sequence[Mapping[str, Any]]],
) -> list[CurvePoint]:
    """Merge spin candidates by taking the lowest-energy point at each distance."""
    best_by_distance: dict[float, CurvePoint] = {}
    for candidate, points in candidate_points.items():
        spin_candidate = str(candidate)
        for point in finite_energy_points(points):
            key = distance_key(point_distance(point))
            energy = point_energy(point)
            current_best = best_by_distance.get(key)
            if current_best is None or energy < point_energy(current_best):
                best_by_distance[key] = {**point, "spin_candidate": spin_candidate}
    return [best_by_distance[key] for key in sorted(best_by_distance)]


def replace_isolated_spin_branch_drops(
    merged_points: Sequence[Mapping[str, Any]],
    candidate_points: Mapping[SpinCandidate, Sequence[Mapping[str, Any]]],
    *,
    min_drop_ev: float = 3.0,
) -> tuple[list[CurvePoint], list[CurvePostprocessEdit]]:
    """Replace isolated severe downward spin-branch drops with the neighbor spin.

    The adiabatic minimum over many spin candidates can contain one-point SCF branch
    failures: an interior point drops by many eV to a different spin candidate, while
    the curve immediately recovers. These isolated downward spikes are not physical
    smooth PEC behavior. If either neighboring spin branch exists at the same distance
    and is closer to the local linear trend, this function replaces the dropped point
    by the smoothest neighboring-branch point.
    """
    processed = [dict(point) for point in merged_points]
    edits: list[CurvePostprocessEdit] = []
    candidate_lookup = _candidate_points_by_spin_and_distance(candidate_points)

    for point_idx in range(1, len(processed) - 1):
        left_point = processed[point_idx - 1]
        point = processed[point_idx]
        right_point = processed[point_idx + 1]

        left_spin = str(left_point.get("spin_candidate", ""))
        point_spin = str(point.get("spin_candidate", ""))
        right_spin = str(right_point.get("spin_candidate", ""))
        if not left_spin or not right_spin or point_spin in {left_spin, right_spin}:
            continue

        left_energy = point_energy(left_point)
        point_energy_ev = point_energy(point)
        right_energy = point_energy(right_point)
        expected_energy = _linear_expected_energy(left_point, point, right_point)
        if point_energy_ev > min(left_energy, right_energy) - min_drop_ev:
            continue
        if point_energy_ev > expected_energy - min_drop_ev:
            continue

        replacement_spin, replacement = _best_neighbor_spin_replacement(
            candidate_lookup,
            point,
            expected_energy,
            left_spin,
            right_spin,
        )
        if replacement is None:
            continue
        replacement_energy = point_energy(replacement)
        original_residual = abs(point_energy_ev - expected_energy)
        replacement_residual = abs(replacement_energy - expected_energy)
        if replacement_energy < point_energy_ev + min_drop_ev:
            continue
        if replacement_residual >= original_residual:
            continue

        processed[point_idx] = {
            **replacement,
            "spin_candidate": replacement_spin,
            "dropped_spin_candidate": point_spin,
            "postprocess_reason": "isolated_spin_branch_energy_drop",
        }
        edits.append(
            CurvePostprocessEdit(
                kind="isolated_spin_branch_drop",
                distance=point_distance(point),
                original_spin=point_spin,
                replacement_spin=replacement_spin,
                original_energy=point_energy_ev,
                replacement_energy=replacement_energy,
                expected_energy=expected_energy,
            )
        )

    return processed, edits


def remove_isolated_energy_bumps(
    merged_points: Sequence[Mapping[str, Any]],
    *,
    min_bump_ev: float = 0.1,
) -> tuple[list[CurvePoint], list[CurvePostprocessEdit]]:
    """Remove isolated upward SCF glitches from otherwise smooth PEC branches."""
    processed = [dict(point) for point in merged_points]
    edits: list[CurvePostprocessEdit] = []
    remove_indices: set[int] = set()
    for point_idx in range(1, len(processed) - 1):
        left_point = processed[point_idx - 1]
        point = processed[point_idx]
        right_point = processed[point_idx + 1]
        point_energy_ev = point_energy(point)
        expected_energy = _linear_expected_energy(left_point, point, right_point)
        if (
            point_energy_ev
            < max(point_energy(left_point), point_energy(right_point)) + min_bump_ev
        ):
            continue
        if point_energy_ev < expected_energy + min_bump_ev:
            continue
        remove_indices.add(point_idx)
        edits.append(
            CurvePostprocessEdit(
                kind="isolated_energy_bump",
                distance=point_distance(point),
                original_spin=str(point.get("spin_candidate", "")),
                replacement_spin=None,
                original_energy=point_energy_ev,
                replacement_energy=None,
                expected_energy=expected_energy,
            )
        )
    return [
        point
        for point_idx, point in enumerate(processed)
        if point_idx not in remove_indices
    ], edits


def _best_neighbor_spin_replacement(
    candidate_lookup: Mapping[tuple[str, float], CurvePoint],
    point: Mapping[str, Any],
    expected_energy: float,
    left_spin: str,
    right_spin: str,
) -> tuple[str, CurvePoint | None]:
    """Return the neighbor spin point closest to the local expected energy."""
    middle_key = distance_key(point_distance(point))
    candidates = [
        (spin, candidate_lookup[(spin, middle_key)])
        for spin in dict.fromkeys((left_spin, right_spin))
        if (spin, middle_key) in candidate_lookup
    ]
    if not candidates:
        return "", None
    return min(
        candidates, key=lambda item: abs(point_energy(item[1]) - expected_energy)
    )


def merge_postprocessed_min_energy_curve(
    candidate_points: Mapping[SpinCandidate, Sequence[Mapping[str, Any]]],
    *,
    min_drop_ev: float = 3.0,
) -> tuple[list[CurvePoint], list[CurvePostprocessEdit]]:
    """Merge spin candidates and remove isolated SCF artifacts."""
    # drop variationally collapsed points per candidate first: they'd always win the
    # per-distance min-merge and no neighbor-branch replacement could outbid them
    cleaned: dict[SpinCandidate, list[CurvePoint]] = {}
    collapse_edits: list[CurvePostprocessEdit] = []
    for candidate, points in candidate_points.items():
        kept, edits = drop_collapsed_scf_points(points, spin_candidate=str(candidate))
        cleaned[candidate] = kept
        collapse_edits += edits

    merged = merge_min_energy_curve(cleaned)
    merged, drop_edits = replace_isolated_spin_branch_drops(
        merged, cleaned, min_drop_ev=min_drop_ev
    )
    merged, bump_edits = remove_isolated_energy_bumps(merged)
    return merged, collapse_edits + drop_edits + bump_edits


def _candidate_points_by_spin_and_distance(
    candidate_points: Mapping[SpinCandidate, Sequence[Mapping[str, Any]]],
) -> dict[tuple[str, float], CurvePoint]:
    """Return finite candidate points keyed by spin label and rounded distance."""
    lookup: dict[tuple[str, float], CurvePoint] = {}
    for candidate, points in candidate_points.items():
        spin_candidate = str(candidate)
        for point in finite_energy_points(points):
            key = (spin_candidate, distance_key(point_distance(point)))
            current = lookup.get(key)
            if current is None or point_energy(point) < point_energy(current):
                lookup[key] = {**point, "spin_candidate": spin_candidate}
    return lookup


def _linear_expected_energy(
    left_point: Mapping[str, Any],
    point: Mapping[str, Any],
    right_point: Mapping[str, Any],
) -> float:
    """Interpolate the local expected energy from immediate neighbors."""
    left_distance = point_distance(left_point)
    middle_distance = point_distance(point)
    right_distance = point_distance(right_point)
    left_energy = point_energy(left_point)
    right_energy = point_energy(right_point)
    if left_distance == right_distance:
        return (left_energy + right_energy) / 2
    fraction = (middle_distance - left_distance) / (right_distance - left_distance)
    return left_energy + fraction * (right_energy - left_energy)
