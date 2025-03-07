# Matbench Discovery Model Submission

## Description

<!-- Please provide a brief description of your model -->

## Checklist

Please check the following items before submitting your PR:

- [ ] I have created a new folder `models/<model_name>` for my submission.
- [ ] I have added the model name to `enums.py`.
- [ ] I have uploaded the required model prediction file to figshare or another cloud storage service (`<yyyy-mm-dd>-<model_name>-preds.csv.gz`).
- [ ] I have uploaded the required optimized geometry file to figshare or another cloud storage service (`<yyyy-mm-dd>-wbm-IS2RE-FIRE.json.gz`).
- [ ] I have uploaded the required phonons file to figshare or another cloud storage service (`<yyyy-mm-dd>-kappa-103-FIRE-dist=0.01-fmax=1e-4-symprec=1e-5.json.gz`).
- [ ] I have included a YAML metadata file (`<model_name>.yml`) with model details and the urls to the figshare files. If not using figshare I have included the urls to the cloud storage service in the description of the PR.
- [ ] I have included the test script (`test_<model_name>_discovery.py`) that generated the dsicovery predictions.

## Additional Information (Optional)

- [ ] I have included a training script (`train_<model_name>.py`) if I trained a model specifically for this benchmark.
- [ ] I have included a README.md with additional details about my model.
