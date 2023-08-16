from glob import glob
from typing import Any

import yaml

from matbench_discovery import ROOT

MODEL_DIRS = glob(f"{ROOT}/models/*/")
MODEL_METADATA: dict[str, dict[str, Any]] = {}

for model_dir in MODEL_DIRS:
    [md_file] = md_files = glob(f"{model_dir}metadata*.yml")
    assert len(md_files) == 1, f"expected 1 metadata file, got {md_files=}"
    if md_file.endswith("aborted.yml"):
        continue
    # make sure all required keys are non-empty
    with open(md_file) as yml_file:
        models = yaml.full_load(yml_file)

    # some metadata files contain a single model, some have multiple
    if not isinstance(models, list):
        models = [models]
    for model in models:
        model["model_dir"] = model_dir
        MODEL_METADATA[model["model_name"]] = model
