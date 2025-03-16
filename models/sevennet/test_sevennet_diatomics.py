"""Generate diatomic potential energy curves for a ML force field."""

# %%
from __future__ import annotations

import gzip
import json
import os
import warnings
from collections import defaultdict
from glob import glob
from typing import Any

import numpy as np
import torch
from sevenn.calculator import SevenNetCalculator

from matbench_discovery import ROOT, today
from matbench_discovery.diatomics import calc_diatomic_curve, homo_nuc

warnings.filterwarnings("ignore", category=DeprecationWarning, module="spglib")


# %% editable config
model_name = "sevennet"
model_variant = "sevennet-mf-ompa"  # choose 7net model variant to eval
device = "cuda" if torch.cuda.is_available() else "cpu"

calculator_kwargs = {
    "sevennet-0": {"model": "7net-0"},
    "sevennet-l3i5": {"model": "7net-l3i5"},
    "sevennet-mf-ompa": {"model": "7net-mf-ompa", "modal": "mpa"},
}["model_variant"]
calculator_kwargs["device"] = device

dtype = "float32"
calculator = SevenNetCalculator(**calculator_kwargs)

json_path = f"{ROOT}/models/{model_name}/{model_variant}/{today}-diatomics.json.gz"
existing_paths = glob(json_path.replace(today, "*-*-*"))
if existing_paths:
    print(f"Skipping {model_name}/{model_variant}\n{existing_paths=}")
    raise SystemExit

# distances in Angstrom
distances = np.linspace(0.1, 6, 100)

# only consider elements up to atomic number 10 (Ne) for now
atomic_nums = range(1, 93)

# generate list of homonuclear pairs [(1, 1), (2, 2), ...]
homo_nuclear_pairs = [(z, z) for z in atomic_nums]

# results structure: model_name->homo_nuc->distances: list[float], energies: list[float]
results: defaultdict[str, dict[str, Any]] = defaultdict(
    lambda: {homo_nuc: {}, "distances": distances}
)

print(f"\nPredicting diatomic curves for {model_name}/{model_variant}")
calc_diatomic_curve(
    pairs=homo_nuclear_pairs,
    calculator=calculator,
    model_name=model_variant,
    distances=distances,
    results=results[model_variant][homo_nuc],
)

out_dir = os.path.dirname(json_path)
os.makedirs(out_dir, exist_ok=True)
print(f"Saving results to {json_path}")

with gzip.open(json_path, mode="wt") as file:
    json.dump(results, file, indent=2, default=lambda x: x.tolist())
