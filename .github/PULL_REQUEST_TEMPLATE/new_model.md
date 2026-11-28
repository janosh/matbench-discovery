## Description

<!-- Please provide a brief description of your model -->

## Checklist

Open PRs as `draft` and mark them `ready to review` only after checking off these items:

- [ ] I created a new folder and YAML metadata file `models/<family>/<model_key>.yml`, where `family` is the lowercase model-family folder and the YAML includes author details, training set names, important hyperparameters, etc.
- [ ] I uploaded all applicable prediction files to Figshare or another cloud storage service:
  - [ ] WBM discovery predictions (`<yyyy-mm-dd>-discovery.csv.gz`).
  - [ ] Model-relaxed structures in [JSON Lines format](https://jsonlines.org) (`<yyyy-mm-dd>-geo-opt.jsonl.gz`).
  - [ ] Normalized 103-material phonon predictions (`<yyyy-mm-dd>-phonons-kappa-103.json.gz`), with any force sets in a separate artifact.
  - [ ] Molecular-dynamics predictions (`<yyyy-mm-dd>-md-metrics.csv.gz`).
  - [ ] Diatomic predictions (`<yyyy-mm-dd>-diatomics.json.gz`).
- [ ] I recorded each file's direct-download URL in the matching `url` key of `models/<family>/<model_key>.yml`.
- [ ] I followed the [shared-runner requirements](https://github.com/janosh/matbench-discovery/blob/main/contributing.md) for calculator registration, kappa settings, and smoke tests.
- [ ] I have run `uv run --with-editable . scripts/ingest_model.py <model_key>` as described in the [contributing guide](https://github.com/janosh/matbench-discovery/blob/main/contributing.md) to check metadata and generate the plots needed for submission.
- [ ] I have installed/run the pre-commit hooks (`prek install` or `uvx prek`) and all checks pass.

Marking your PR `ready to review` triggers an automated code review; please address the issues it raises. Minimizing the final human review lets maintainers merge your submission much faster!

## Additional Information (Optional)

- [ ] I included a training script (`train_<model_key>.py` or `train_<family>.py`) if I trained a model specifically for this benchmark.
- [ ] I included a `readme.md` with additional details about my model.
