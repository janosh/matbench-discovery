# /// script
# requires-python = ">=3.14"
# dependencies = ["numpy>=2.3", "ase>=3.26", "upet==0.3.1", "matbench-discovery"]
# [tool.uv.sources]
# matbench-discovery = { path = "../../", editable = true }
# ///
"""Reproduce the R24E elemental energy references from the packaged development set."""

from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path

import numpy as np

DIRECTORY = Path(__file__).resolve().parent / "pet-omat-xs-r24e"


def fit(
    rows: list[dict], orientation_energies: list[float]
) -> tuple[float, np.ndarray]:
    """Fit the per-atom intercept and a ridge-regularized per-element reference."""
    counts = Counter(row["formula"] for row in rows)
    weights = np.array([1 / (len(counts) * counts[row["formula"]]) for row in rows])
    composition = np.asarray(
        [
            np.bincount(row["atomic_numbers"], minlength=119)
            / len(row["atomic_numbers"])
            for row in rows
        ]
    )
    x = np.column_stack((np.ones(len(rows)), composition))
    target = np.asarray(
        [
            row["dft_energy"] / len(row["atomic_numbers"])
            - prediction / len(row["atomic_numbers"])
            for row, prediction in zip(rows, orientation_energies, strict=True)
        ]
    )
    theta = np.zeros(x.shape[1])
    theta[0] = np.median(target)
    penalty = np.diag([0.0] + [0.0005] * 119)
    for _ in range(6):
        w = weights / np.maximum(1.0, abs(target - x @ theta) / 0.1)
        theta = np.linalg.solve((x.T * w) @ x + penalty, (x.T * w) @ target)
    return float(theta[0]), theta[1:]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Optional fitted reference JSON")
    parser.add_argument(
        "--recompute-parent",
        action="store_true",
        help="Recompute all 24 orientation energies",
    )
    parser.add_argument(
        "--checkpoint", help="Optional local PET-OMat-XS checkpoint override"
    )
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cpu")
    args = parser.parse_args()
    with gzip.open(
        DIRECTORY / "training-energies.jsonl.gz", "rt", encoding="utf-8"
    ) as f:
        rows = [json.loads(line) for line in f]
    frozen = json.loads((DIRECTORY / "frozen-model.json").read_text(encoding="utf-8"))
    if [r["id"] for r in rows] != frozen["fit_ids"] or len(
        set(frozen["fit_ids"])
    ) != 871:
        raise ValueError("Dataset changed relative to the released fit")
    predictions = [r["r24_energy"] for r in rows]
    if args.recompute_parent:
        from ase import Atoms

        from matbench_discovery.calculators import load_calculator
        from matbench_discovery.pet_r24e import R24E

        registered = load_calculator(
            "pet-omat-xs-r24e",
            device=args.device,
            dtype="float32",
            checkpoint=args.checkpoint,
        )
        # Raw R24 energies, before learning the shift. The calculator is intentionally
        # sequential to use the reference implementation's ordering on all systems.
        parent = R24E(registered.base, intercept=0.0, offsets=np.zeros(119))
        predictions = []
        for row in rows:
            atoms = Atoms(
                numbers=row["atomic_numbers"],
                cell=row["cell"],
                scaled_positions=row["scaled_positions"],
                pbc=True,
                calculator=parent,
            )
            predictions.append(atoms.get_potential_energy())
    intercept, offsets = fit(rows, predictions)
    # Native float32 PET evaluation can differ slightly with backend/reduction order.
    tolerance = 3e-5 if args.recompute_parent else 1e-10
    np.testing.assert_allclose(
        [intercept, *offsets],
        [frozen["intercept"], *frozen["offsets"]],
        atol=tolerance,
        rtol=tolerance,
    )
    if args.output is not None:
        args.output.write_text(
            json.dumps(dict(intercept=intercept, offsets=offsets.tolist()), indent=2)
            + "\n"
        )
    print(f"Reproduced frozen references from {len(rows)} development structures.")


if __name__ == "__main__":
    main()
