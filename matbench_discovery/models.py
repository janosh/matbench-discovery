"""Initializes global variable MODEL_METADATA."""

from glob import glob
from typing import Any

import yaml

from matbench_discovery import ROOT

# ignore underscore-prefixed directories for WIP model submissions
MODEL_DIRS = glob(f"{ROOT}/models/[!_]*/")
MODEL_METADATA: dict[str, dict[str, Any]] = {}

for model_dir in MODEL_DIRS:
    md_files = glob(f"{model_dir}*.yml")
    if not 1 <= len(md_files) <= 2:
        raise RuntimeError(f"expected 1 metadata file, got {md_files=} in {model_dir=}")
    for md_file in md_files:
        if md_file.endswith("aborted.yml"):
            continue
        # make sure all required keys are non-empty
        with open(md_file) as yml_file:
            model_data = yaml.full_load(yml_file)

        model_data["model_dir"] = model_dir
        MODEL_METADATA[model_data["model_name"]] = model_data
