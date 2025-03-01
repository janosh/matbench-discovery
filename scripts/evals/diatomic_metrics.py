"""Calculate diatomic curve metrics for all models and write them to YAML files."""

import gzip
import json
import os

from matbench_discovery import ROOT
from matbench_discovery.enums import Model
from matbench_discovery.metrics import diatomics
from matbench_discovery.metrics.diatomics import DiatomicCurves

# Loop over all models
for model in Model:
    if not os.path.isfile(model.yaml_path):
        continue

    diatomics_metrics = model.metrics.get("diatomics")
    if not isinstance(diatomics_metrics, dict):
        continue

    pred_file = diatomics_metrics.get("pred_file")
    if not isinstance(pred_file, str):
        continue

    abs_path = f"{ROOT}/{pred_file}"
    if not os.path.isfile(abs_path):
        print(f"Prediction file {pred_file} not found for {model.name}")
        continue

    # Load predicted curves
    with gzip.open(abs_path, mode="rb") as file:
        pred_data = json.load(file) or {}

    pred_curves = DiatomicCurves.from_dict(pred_data)

    # Calculate metrics (without reference data)
    metrics = diatomics.calc_diatomic_metrics(ref_curves=None, pred_curves=pred_curves)

    # Write metrics to YAML
    diatomics.write_metrics_to_yaml(model, metrics)
