"""Build the bundled DFT diatomic reference from spin-candidate VASP curves."""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
from dataclasses import asdict
from typing import Any

from matbench_discovery import ROOT
from matbench_discovery.metrics.diatomics.reference import (
    CurvePostprocessEdit,
    count_dissociation_tail_jumps,
    count_magmom_discontinuities,
    merge_min_energy_curve,
    merge_postprocessed_min_energy_curve,
    point_distance,
    point_energy,
)

FUNC_LABEL = {"pbe": "PBE", "r2scan": "r2SCAN"}
XC_ORDER = ("pbe", "r2scan")
MIN_COMPLETE_POINTS = 45


def has_finite_forces(forces: object) -> bool:
    """Return whether forces are two finite Cartesian 3-vectors."""
    if not isinstance(forces, list) or len(forces) != 2:
        return False
    return all(
        isinstance(atom, list)
        and len(atom) == 3
        and all(
            isinstance(component, int | float) and math.isfinite(component)
            for component in atom
        )
        for atom in forces
    )


def spin_suffix(candidate: int | str) -> str:
    """Return the output directory suffix for one spin candidate token."""
    return "n0afm" if candidate == "afm" else f"n{candidate}"


def load_candidate_points(
    src_dir: str, symbol: str, xc: str, candidate: int | str
) -> list[dict[str, Any]]:
    """Load one candidate curve's finite-energy and finite-force points."""
    path = f"{src_dir}/{symbol}_{xc}_{spin_suffix(candidate)}/curve.json"
    if not os.path.isfile(path):
        return []
    with open(path) as file:
        points = json.load(file).get("points", [])

    finite_points: list[dict[str, Any]] = []
    for point in points:
        try:
            point_energy(point)
            point_distance(point)
            if not has_finite_forces(point.get("forces")):
                continue
        except (KeyError, TypeError, ValueError):
            continue
        finite_points.append(point)
    return finite_points


def serializable_curve(points: list[dict[str, Any]]) -> dict[str, list]:
    """Convert merged curve points to the bundled reference JSON schema.

    Besides distances/energies/forces, each point carries its site-projected magnetic
    moments (LORBIT=11, eV-scale spin-state hops show up as magmom discontinuities)
    and the winning spin candidate ("<NUPDOWN>" or "afm") so spin-aware MLIPs can be
    compared against the reference spin state directly (suggested by Andrew Rosen).
    """
    distances: list[float] = []
    energies: list[float] = []
    forces: list[list[list[float]]] = []
    magmoms: list[list[float] | None] = []
    spin_candidates: list[str | None] = []
    for point in points:
        distances.append(round(point["distance"], 6))
        energies.append(round(point["energies"][-1], 6))
        forces.append([[round(comp, 6) for comp in atom] for atom in point["forces"]])
        point_magmoms = point.get("magmoms")
        magmoms.append(
            [round(moment, 3) for moment in point_magmoms] if point_magmoms else None
        )
        spin_candidate = point.get("spin_candidate")
        spin_candidates.append(None if spin_candidate is None else str(spin_candidate))
    return {
        "distances": distances,
        "energies": energies,
        "forces": forces,
        "magmoms": magmoms,
        "spin_candidates": spin_candidates,
    }


def build_reference(
    *,
    src_dir: str,
    candidate_map_path: str,
    out_path: str,
    merged_dir: str | None,
    min_drop_ev: float,
    postprocess: bool,
) -> dict[str, dict[str, int]]:
    """Build and write the reference file, returning summary counts."""
    with open(candidate_map_path) as file:
        candidate_map: dict[str, list[int | str]] = json.load(file)

    refs: dict[str, dict[str, dict[str, list]]] = {"PBE": {}, "r2SCAN": {}}
    summary = {
        "merged": {},
        "skipped": {},
        "short_candidate_pairs": {},
        "tail_jump_pairs": {},
        "magmom_jump_pairs": {},
        "postprocess_edits": {},
    }
    quality_rows: list[dict[str, object]] = []
    for symbol, candidates in candidate_map.items():
        for xc in XC_ORDER:
            candidate_points = {
                candidate: load_candidate_points(src_dir, symbol, xc, candidate)
                for candidate in candidates
            }
            short_candidates = {
                candidate: len(points)
                for candidate, points in candidate_points.items()
                if len(points) < MIN_COMPLETE_POINTS
            }
            if postprocess:
                merged_points, replacements = merge_postprocessed_min_energy_curve(
                    candidate_points, min_drop_ev=min_drop_ev
                )
            else:
                merged_points = merge_min_energy_curve(candidate_points)
                replacements = []
            label = FUNC_LABEL[xc]
            if len(merged_points) < 2:
                # surface dropped elements in the summary and quality report: a pair
                # silently vanishing from the bundled reference between rebuilds
                # (e.g. all candidate curves missing/corrupt) should be visible
                summary["skipped"][label] = summary["skipped"].get(label, 0) + 1
                quality_rows.append(
                    {
                        "symbol": symbol,
                        "xc": xc,
                        "merged_points": len(merged_points),
                        "short_candidates": short_candidates,
                        "skipped": True,
                    }
                )
                continue

            refs[label][f"{symbol}-{symbol}"] = serializable_curve(merged_points)
            summary["merged"][label] = summary["merged"].get(label, 0) + 1
            summary["postprocess_edits"][label] = summary["postprocess_edits"].get(
                label, 0
            ) + len(replacements)
            tail_jumps = count_dissociation_tail_jumps(
                [point_energy(point) for point in merged_points]
            )
            magmom_jumps = count_magmom_discontinuities(
                [point.get("magmoms") for point in merged_points]
            )
            if tail_jumps:
                summary["tail_jump_pairs"][label] = (
                    summary["tail_jump_pairs"].get(label, 0) + 1
                )
            if magmom_jumps:
                summary["magmom_jump_pairs"][label] = (
                    summary["magmom_jump_pairs"].get(label, 0) + 1
                )
            if short_candidates:
                summary["short_candidate_pairs"][label] = (
                    summary["short_candidate_pairs"].get(label, 0) + 1
                )
            if short_candidates or tail_jumps or magmom_jumps:
                quality_rows.append(
                    {
                        "symbol": symbol,
                        "xc": xc,
                        "merged_points": len(merged_points),
                        "short_candidates": short_candidates,
                        "tail_jumps": tail_jumps,
                        "magmom_jumps": magmom_jumps,
                    }
                )
            if merged_dir:
                write_merged_curve(merged_dir, symbol, xc, merged_points, replacements)

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with gzip.open(out_path, "wt", encoding="utf-8") as file:
        json.dump(refs, file, separators=(",", ":"))
    if merged_dir:
        os.makedirs(merged_dir, exist_ok=True)
        with open(f"{merged_dir}/reference-quality.json", "w") as file:
            json.dump(quality_rows, file, indent=2)
    return summary


def write_merged_curve(
    merged_dir: str,
    symbol: str,
    xc: str,
    merged_points: list[dict[str, Any]],
    replacements: list[CurvePostprocessEdit],
) -> None:
    """Write one merged diagnostic curve with postprocess metadata."""
    out_dir = f"{merged_dir}/{symbol}_{xc}"
    os.makedirs(out_dir, exist_ok=True)
    with open(f"{out_dir}/curve.json", "w") as file:
        postprocess_data = {
            "postprocess_edits": [asdict(edit) for edit in replacements]
        }
        json_data = {
            "element": symbol,
            "functional": xc,
            "points": merged_points,
            "postprocess": postprocess_data,
        }
        json.dump(json_data, file)


def main() -> None:
    """Parse CLI arguments and build the DFT diatomic reference."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--src-dir", default=f"{ROOT}/tmp/diatomics-candidates")
    parser.add_argument("--candidate-map", default=f"{ROOT}/tmp/adiabatic_cands.json")
    parser.add_argument(
        "--out-path", default=f"{ROOT}/site/src/lib/diatomics-dft.json.gz"
    )
    parser.add_argument("--merged-dir", default=f"{ROOT}/tmp/diatomics-merged")
    parser.add_argument("--min-drop-ev", type=float, default=3.0)
    parser.add_argument("--no-postprocess", action="store_true")
    args = parser.parse_args()

    summary = build_reference(
        src_dir=args.src_dir,
        candidate_map_path=args.candidate_map,
        out_path=args.out_path,
        merged_dir=args.merged_dir,
        min_drop_ev=args.min_drop_ev,
        postprocess=not args.no_postprocess,
    )
    size_kb = os.path.getsize(args.out_path) / 1024
    print(f"wrote {args.out_path} ({size_kb:.0f} KB)")
    print(summary)


if __name__ == "__main__":
    main()
