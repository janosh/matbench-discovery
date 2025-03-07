## Description

<!-- Please provide a brief description of your model -->

## Checklist

Please check the following items before submitting your PR:

- [ ] I have created a new folder `models/<arch_name>/<model_variant>` for my submission. `arch_name` is the name of the architecture and `model_variant` includes things like training set names and/or important hyperparameters.
- [ ] I have added the my new model as a new attribute on the [`Model.<arch_name>` enum](https://github.com/janosh/matbench-discovery/blob/57d0d0c8a14cd317/matbench_discovery/enums.py#L274) in `enums.py`.
- [ ] I have uploaded the energy/force/stress model prediction file for the WBM test set to Figshare or another cloud storage service (`<yyyy-mm-dd>-<arch_name>-preds.csv.gz`).
- [ ] I have uploaded the model-relaxed structures file to Figshare or another cloud storage service (`<yyyy-mm-dd>-wbm-IS2RE-FIRE.json.gz`).
- [ ] I have uploaded the phonon predictions to Figshare or another cloud storage service (`<yyyy-mm-dd>-kappa-103-FIRE-<values-of-dist|fmax|symprec>.gz`).
- [ ] I have included a YAML metadata file (`models/<arch_name>/<model_variant>.yml`) with model details and the urls to the Figshare files. If not using Figshare I have included the urls to the cloud storage service in the description of the PR.
- [ ] I have included the test script (`test_<arch_name>_<task>.py` for `task` in `discovery`, `phonons`) that generated the prediction files.

## Additional Information (Optional)

- [ ] I have included a training script (`train_<arch_name>.py`) if I trained a model specifically for this benchmark.
- [ ] I have included a `readme.md` with additional details about my model.
