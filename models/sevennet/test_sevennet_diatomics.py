"""Generate diatomic potential energy curves for a ML force field."""

# %%
from __future__ import annotations

import gzip
import json
import os
import warnings
from collections import defaultdict
from glob import glob
from pathlib import Path
from typing import Any

import numpy as np
import torch
from sevenn.calculator import SevenNetCalculator

from matbench_discovery import ROOT, today
from matbench_discovery.diatomics import calc_diatomic_curve, homo_nuc

warnings.filterwarnings("ignore", category=DeprecationWarning, module="spglib")


# %% editable config
model_name = "sevennet"
model_variant = "sevennet-omni-i12"  # choose 7net model variant to eval
device = "cuda" if torch.cuda.is_available() else "cpu"

calc_kwargs: dict[str, Any] = {
    "sevennet-0": {"model": "7net-0"},
    "sevennet-l3i5": {"model": "7net-l3i5"},
    "sevennet-mf-ompa": {"model": "7net-mf-ompa", "modal": "mpa"},
    "sevennet-omni-i12": {"model": "7net-omni-i12", "modal": "mpa"},
}[model_variant]
calc_kwargs["device"] = device

# Will be removed after integrating model checkpoint download into sevenn package
checkpoint_urls = {
    "sevennet-omni-i12": "https://figshare.com/ndownloader/files/60977863",
}
if model_variant in checkpoint_urls:
    cache_dir = Path.home() / ".cache" / "sevennet"
    cache_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_path = cache_dir / f"checkpoint_{model_variant.replace('-', '_')}.pth"

    if not checkpoint_path.exists():
        print(f"Downloading {model_variant} checkpoint to {checkpoint_path}...")
        import requests

        response = requests.get(checkpoint_urls[model_variant], stream=True, timeout=30)
        response.raise_for_status()
        with open(checkpoint_path, "wb") as f:
            f.writelines(response.iter_content(chunk_size=8192))
        print("Download complete.")
    else:
        print(f"Using cached checkpoint: {checkpoint_path}")
    calc_kwargs["model"] = str(checkpoint_path)

dtype = "float32"
calculator = SevenNetCalculator(**calc_kwargs)

json_path = f"{ROOT}/models/{model_name}/{model_variant}/{today}-diatomics.json.gz"
existing_paths = glob(json_path.replace(today, "*-*-*"))
if existing_paths:
    print(f"Skipping {model_name}/{model_variant}\n{existing_paths=}")
    raise SystemExit


# log-spaced atomic separations going inwards from 10 to 0.1 Angstrom
distances = np.logspace(1, -1, 40)

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
