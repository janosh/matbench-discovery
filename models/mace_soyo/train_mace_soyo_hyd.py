# /// script
# requires-python = ">=3.14"
# dependencies = ["numpy==2.5.4"]
# ///
"""Refit the frozen 16-parameter correction from independent heteronuclear data.

Run from the repository root after the HYD package has been added to it:
  python models/mace_soyo/train_mace_soyo_hyd.py --output refit-candidate.json
This repeats the selected ridge fit, not benchmark-driven parameter selection.
"""

import argparse
import hashlib
import importlib
import json
import sys
import zipfile
from pathlib import Path

import numpy as np


def main(output: Path) -> None:
    """Repeat the selected independent-data ridge fit and verify frozen parameters."""
    folder = Path(__file__).with_name("mace-soyo-hyd-radial-v8")
    archive = folder / "frozen-source.zip"
    source_hash = "f29f5d21250764e3b5e4128521ca118a4e988f54db51138ea661369c5061ab2c"
    if hashlib.sha256(archive.read_bytes()).hexdigest() != source_hash:
        raise ValueError("The original scientific source archive failed its checksum")
    sys.path.insert(0, str(archive))
    transfer = importlib.import_module("fit_transfer_candidate")
    fit = importlib.import_module("screen_joint_balance").fit
    names = ["initial-fit.jsonl", "transfer-fit.jsonl"]
    hashes = [
        "7aedca2539e1421a1a4a021cf0d7588261acac9c010de3d84ef0d1c1ab14f0a7",
        "b21fc8ca4bb51bbb32c551e76ffc9137880ca3b45472eebf535070b883ab0e6e",
    ]
    with zipfile.ZipFile(archive) as data:
        paths = [zipfile.Path(data, "training/" + name) for name in names]
        for path, expected in zip(paths, hashes, strict=True):
            if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                raise ValueError(f"Training file checksum mismatch: {path}")
        groups = {
            pair: rows
            for pair, rows in transfer.collect(paths).items()
            if len(rows) >= 5
        }
    if len(groups) != 24:
        raise ValueError("Expected exactly the 24 eligible independent training groups")
    designs = {
        pair: transfer.design(pair, rows, "radial") for pair, rows in groups.items()
    }
    coefficients = fit(designs, groups, penalty=0.85, energy_weight=16.0)
    frozen = (
        Path(__file__).resolve().parents[2]
        / "matbench_discovery/mace_soyo_hyd/frozen-candidate.json"
    )
    expected = np.array(json.loads(frozen.read_text())["coefficients"])
    if not np.allclose(coefficients, expected, atol=1e-10, rtol=1e-10):
        raise ValueError(
            "Refitted parameters failed the predeclared frozen-model tolerance"
        )
    result = {
        "family": "radial",
        "coefficients": coefficients.tolist(),
        "ridge_penalty": 0.85,
        "energy_weight": 16.0,
        "training_molecule_groups": len(groups),
        "max_coefficient_absolute_difference": float(
            np.max(np.abs(coefficients - expected))
        ),
        "training_sha256": dict(zip(names, hashes, strict=True)),
        "official_labels_used": False,
    }
    with output.open("x") as stream:
        stream.write(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                "groups": len(groups),
                "max_coefficient_difference": result[
                    "max_coefficient_absolute_difference"
                ],
            }
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    main(parser.parse_args().output)
