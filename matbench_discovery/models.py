"""Initializes global variable MODEL_METADATA."""

from glob import glob
from typing import Any

import yaml

from matbench_discovery import ROOT
from matbench_discovery.enums import ModelType, Open

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
        # make sure all required keys are non-empty
        with open(metadata_file) as yml_file:
            model_data = yaml.full_load(yml_file)

        # skip models that aren't completed
        if model_data.get("status", "complete") != "complete":
            continue

        model_data["model_dir"] = model_dir
        # YAML file name serves as model key and must match Model enum attribute
        model_data["model_key"] = metadata_file.split("/")[-1].split(".yml")[0]
        MODEL_METADATA[model_data["model_name"]] = model_data

        try:
            ModelType(model_data.get("model_type"))  # check if model_type is valid
        except ValueError as exc:
            exc.add_note(f"{metadata_file=}")
            raise


def model_is_compliant(metadata: dict[str, str | list[str]]) -> bool:
    """Check if model complies with benchmark restrictions on allowed training sets and
    open model code and weights.

    Args:
        metadata: model metadata dictionary

    Returns:
        bool: True if model is compliant
    """
    openness = metadata.get("openness", Open.OSOD)
    if openness != Open.OSOD:
        return False
    training_sets = metadata.get("training_set")

    if not isinstance(training_sets, list):
        model_name = metadata.get("model_name")
        raise TypeError(
            f"{model_name}: expected list of training sets, got {training_sets=}"
        )

    return set(training_sets) <= {"MP 2022", "MPtrj", "MPF", "MP Graphs"}
