## ALIGNN formation energy predictions on WBM test set

ALIGNN is trained using L1 loss and 1000 epochs. The model that performs best on the validation set is saved and used for predictions (requires minor adaptation of the ALIGNN source). The modifications to the ALIGNN source code are provided as patch `alignn-2023.01.10.patch`, which was applied to ALIGNN version `2023.01.10`. In addition, all Python requirements are given as `requirements.txt`.

To reproduce the `alignn` package state used for this submission, run

```bash
pip install alignn==2023.01.10
alignn_dir=$(python -c "import alignn; print(alignn.__path__[0])")
cd $alignn_dir
git apply /path/to/alignn-2023.01.10.patch
```

Replace `/path/to/` with the actual path to the patch file.

The directory contains the following files, which must be executed in the given order to reproduce the results:

1. `train_data.py`: Export Matbench Discovery training data to ALIGNN compatible format. This script outputs training data in the directory `data_train`. In addition, a small test data set is set apart and stored in the directory `data_test`
1. `train_alignn.py`: Train an ALIGNN model on previously exported data. The resulting model is stored in the directory `data-train-result`
1. `test_data.py`: Export WBM test data in ALIGNN-compatible format. The data is stored in the directory `data-test-wbm`
1. `test_alignn.py`: Test a trained ALIGNN model on the WBM data. Predictions are stored in the file `test_alignn_result.json`
