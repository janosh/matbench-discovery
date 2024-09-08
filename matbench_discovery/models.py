"""Initializes global variable MODEL_METADATA."""

from glob import glob
from typing import Any

import yaml

from matbench_discovery import ROOT
from matbench_discovery.enums import ModelType

# ignore underscore-prefixed directories for WIP model submissions
MODEL_DIRS = glob(f"{ROOT}/models/[!_]*/")
MODEL_METADATA: dict[str, dict[str, Any]] = {}

for model_dir in MODEL_DIRS:
    metadata_files = glob(f"{model_dir}*.yml")
    if not 1 <= len(metadata_files) <= 2:
        raise RuntimeError(
            f"expected 1 metadata file, got {metadata_files=} in {model_dir=}"
        )
    for metadata_file in metadata_files:
        if metadata_file.endswith("aborted.yml"):
            continue
        # make sure all required keys are non-empty
        with open(metadata_file) as yml_file:
            model_data = yaml.full_load(yml_file)

        model_data["model_dir"] = model_dir
        # YAML file name serves as model key and must match Model enum attribute
        model_data["model_key"] = metadata_file.split("/")[-1].split(".yml")[0]
        MODEL_METADATA[model_data["model_name"]] = model_data

        try:
            ModelType(model_data.get("model_type"))  # check if model_type is valid
        except ValueError as exc:
            exc.add_note(f"{metadata_file=}")
            raise
