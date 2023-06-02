## ALIGNN formation energy predictions on WBM test set

ALIGNN is trained using L1 loss and 1000 epochs. The model that performs best on the validation set is saved and used for predictions (requires minor adaptation of the ALIGNN source). The modifications to the ALIGNN source code are provided as patch `train_alignn_patch.txt`, which was applied to ALIGNN tag `v2023.01.10`. In addition, all python requirements are given as `requirements.txt`.

The directory contains the following files, which must be executed in the given order to reproduce the results:
* `train_data.py`: Export matbench discovery training data to ALIGNN compatible format. This script outputs training data in the directory `data_train`. In addition, a small test data set is set apart and stored in the directory `data_test`
* `train_alignn.py`: Train an ALIGNN model on previously exported data. The resulting model is stored in the directory `data_train_result`
* `test_data.py`: Export WBM test data in ALIGNN compatible format. The data is stored in the directory `data_test_wbm`
* `test_alignn.py`: Test a trained ALIGNN model on the WBM data. Predictions are stored in the file `test_alignn_result.json`
