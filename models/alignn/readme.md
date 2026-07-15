# ALIGNN

Atomistic Line Graph Neural Network for direct formation-energy prediction.

## ALIGNN formation energy predictions on WBM test set

ALIGNN was trained for 1000 epochs using L1 loss. The model that performed best on the validation set was [uploaded to Figshare](https://figshare.com/articles/dataset/Matbench_Discovery_v1_0_0/22715158?file=41233560) as [`2023-06-02-pbenner-best-alignn-model.pth.zip`](https://figshare.com/files/41233560) and used for predictions. This required minor changes to the ALIGNN source code provided in [`alignn/alignn-2023.01.10.patch`](alignn/alignn-2023.01.10.patch).

1. Fix use without test set (see [ALIGNN #104)](https://github.com/usnistgov/alignn/issues/104#issue-1723978225). In this case, we forked a test set, but it might be better to use the entire data, as mentioned above, especially if the test set by chance contains some important outliers.
1. The `Checkpoint` handler in ALIGNN does not define a score name (see [`train.py`](https://github.com/usnistgov/alignn/blob/46334500/alignn/train.py#L851)), so it will just save the last two models during training. With this patch, also the best model in terms of accuracy on the validation set is saved, which is the one used to make predictions. This is important because I used a relatively large `n_early_stopping` in case the validation accuracy shows a double descent (see [Figure 10](https://arxiv.org/pdf/1912.02292.pdf)).

The changes in `alignn/alignn-2023.01.10.patch` were applied to ALIGNN version [`2023.01.10`](https://pypi.org/project/alignn/2023.1.10).

To reproduce the `alignn` package state used for this submission, run

```bash
pip install alignn==2023.01.10
alignn_dir=$(python -c "import alignn; print(alignn.__path__[0])")
cd $alignn_dir
git apply /path/to/models/alignn/alignn/alignn-2023.01.10.patch
```

Replace `/path/to/` with the actual path to the patch file.

The directory retains the training workflow and metadata used for the published
predictions:

1. [`train_alignn.py`](train_alignn.py): Train an ALIGNN model on all 154k MP-computed structure entries. Run it with `uv run models/alignn/train_alignn.py`; its inline metadata installs ALIGNN, W&B, and this repository in an isolated environment. The script uses the original [`alignn-config.json`](alignn-config.json), ported to the current ALIGNN data-loader API. The resulting model checkpoint is saved to the `out_dir` variable in that script; the published checkpoint is linked above.

The legacy one-off discovery executable that generated
`2023-06-03-mp-e-form-alignn-wbm-IS2RE.csv.gz` was retired when
[PR #375](https://github.com/janosh/matbench-discovery/pull/375) centralized
discovery execution in `models/run_discovery.py`. ALIGNN is an archived
direct-prediction submission and does not expose an ASE calculator for that runner.
