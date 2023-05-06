import os
from glob import glob

import yaml

from matbench_discovery import ROOT

MODEL_DIRS = glob(f"{ROOT}/models/*/")


def test_model_dirs_have_metadata() -> None:
    required = (
        "authors",
        "date_added",
        "matbench_discovery_version",
        "model_name",
        "model_version",
        "repo",
    )
    for model_dir in MODEL_DIRS:
        md_file = f"{model_dir}metadata.yml"
        assert os.path.isfile(md_file), f"Missing metadata file: {md_file}"

        # make sure all required keys are non-empty
        with open(md_file) as yml_file:
            models = yaml.full_load(yml_file)

        # some metadata files have a single model, some have multiple
        if not isinstance(models, list):
            models = [models]

        for metadata in models:
            for key in required:
                assert metadata.get(key), f"Empty required {key=} in {md_file}"

            authors, date_added, mbd_version, model_name, model_version, repo = (
                metadata[key] for key in required
            )

            # make sure all keys are valid
            for name in model_name if isinstance(model_name, list) else [model_name]:
                assert (
                    3 < len(name) < 50
                ), f"Invalid {name=} not between 3 and 50 characters"
            assert (
                1 < len(model_version) < 15
            ), f"Invalid {model_version=} not between 1 and 15 characters"
            # TODO increase max allowed version when updating package
            assert (
                1 <= mbd_version <= 1
            ), f"Invalid matbench-discovery version: {mbd_version}"
            assert isinstance(date_added, str), f"Invalid {date_added=} not a string"
            assert isinstance(authors, list)
            assert 1 < len(authors) < 30, f"{len(authors)=} not between 1 and 30"
            assert repo.startswith(
                "https://"
            ), f"Invalid {repo=} not starting with https://"


def test_model_dirs_have_test_scripts() -> None:
    for model_dir in MODEL_DIRS:
        test_scripts = glob(f"{model_dir}*test_*.py")
        test_nbs = glob(f"{model_dir}*test_*.ipynb")
        assert len(test_scripts + test_nbs) > 0, f"Missing test file in {model_dir}"
