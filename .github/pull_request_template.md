## Description

<!-- Please provide a brief description of your model -->

## Checklist

Please open PRs as `draft` and only mark as `ready to review` after checking off the following items:

- [ ] I created a new folder and YAML metadata file `models/<arch_name>/<model_variant>.yml` for my submission. `arch_name` is the name of the architecture and `model_variant.yml` includes things like author details, training set names and important hyperparameters.
- [ ] I added the model to the [`Model` enum](https://github.com/janosh/matbench-discovery/blob/main/matbench_discovery/enums.py).
- [ ] I uploaded the energy/force/stress model prediction file for the WBM test set to Figshare or another cloud storage service (`<yyyy-mm-dd>-<model_variant>-preds.csv.gz`).
- [ ] I uploaded the model-relaxed structures file to Figshare or another cloud storage service in [JSON lines format](https://jsonlines.org) (`<yyyy-mm-dd>-wbm-IS2RE-FIRE.jsonl.gz`). JSON Lines allows fast loading of small numbers of structures with `pandas.read_json(lines=True, nrows=100)` for inspection.
- [ ] I uploaded the normalized 103-material phonon predictions to Figshare or another cloud storage service (`<yyyy-mm-dd>-phonondb-kappa-103.json.gz`), with any force sets in a separate artifact.
- [ ] I have uploaded the diatomic predictions to Figshare or another cloud storage service (`<yyyy-mm-dd>-diatomics.json.gz`).
- [ ] I included the urls to the Figshare files in the YAML metadata file (`models/<arch_name>/<model_variant>.yml`). If not using Figshare I have included the urls to the cloud storage service in the description of the PR.
- [ ] I followed the [shared-runner requirements](https://github.com/janosh/matbench-discovery/blob/main/contributing.md) for calculator registration, kappa settings, and smoke tests.
- [ ] I have run `uv run --with-editable . scripts/ingest_model.py <model_variant>` as described in the [contributing guide](https://github.com/janosh/matbench-discovery/blob/main/contributing.md) to check metadata and generate the plots needed for submission.
- [ ] I have installed/run the pre-commit hooks (`prek install` or `uvx prek`) and ensured all checks are passing.

Marking your PR as `ready to review` will trigger an automated code review, please address any issues raised by the automated review. For the maintainers minimizing the final human review process enables us to merge your submissions much faster!

## Additional Information (Optional)

- [ ] I included a training script (`train_<arch_name>.py`) if I trained a model specifically for this benchmark.
- [ ] I included a `readme.md` with additional details about my model.
