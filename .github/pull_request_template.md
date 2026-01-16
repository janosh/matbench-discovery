## Description

<!-- Please provide a brief description of your model -->

## Checklist

Please check the following items before submitting your PR:

- [ ] I created a new folder and YAML metadata file `models/<arch_name>/<model_variant>.yml` for my submission. `arch_name` is the name of the architecture and `model_variant.yml` includes things like author details, training set names and important hyperparameters.
- [ ] I added the new model as a new attribute [`Model.<model_variant>` enum](https://github.com/janosh/matbench-discovery/blob/57d0d0c8a14cd317/matbench_discovery/enums.py#L274) on the `Model` enum in `enums.py`.
- [ ] I uploaded the energy/force/stress model prediction file for the WBM test set to Figshare or another cloud storage service (`<yyyy-mm-dd>-<model_variant>-preds.csv.gz`).
- [ ] I uploaded the model-relaxed structures file to Figshare or another cloud storage service in [JSON lines format](https://jsonlines.org) (`<yyyy-mm-dd>-wbm-IS2RE-FIRE.jsonl.gz`). JSON Lines allows fast loading of small numbers of structures with `pandas.read_json(lines=True, nrows=100)` for inspection.
- [ ] I uploaded the phonon predictions to Figshare or another cloud storage service (`<yyyy-mm-dd>-kappa-103-FIRE-<values-of-dist|fmax|symprec>.gz`).
- [ ] I have uploaded the diatomic predictions to Figshare or another cloud storage service (`<yyyy-mm-dd>-diatomics.json.gz`).
- [ ] I included the urls to the Figshare files in the YAML metadata file (`models/<arch_name>/<model_variant>.yml`). If not using Figshare I have included the urls to the cloud storage service in the description of the PR.
- [ ] I included the test script (`test_<arch_name>_<task>.py` for `task` in `discovery`, `kappa`, `diatomics`) that generated the prediction files.
- [ ] I have run `just prepare-model-submission <model_variant>` as described in the [contributing guide](https://github.com/janosh/matbench-discovery/contributing.md) to check metadata and generate the plots needed for submission.
- [ ] I have installed/run the pre-commit hooks (`prek install` or `uvx prek`)

## Additional Information (Optional)

- [ ] I included a training script (`train_<arch_name>.py`) if I trained a model specifically for this benchmark.
- [ ] I included a `readme.md` with additional details about my model.
