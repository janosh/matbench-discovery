## MACE formation energy predictions on WBM test set

This submission uses the [`2023-07-14-mace-universal-2-big-128-6.model`](https://figshare.com/ndownloader/files/41565618) checkpoint pre-trained on the [original M3GNet training set](https://figshare.com/articles/dataset/MPF_2021_2_8/19470599).

MACE relaxed each test set structure until the maximum force in the training set dropped below 0.05 eV/Å or 500 optimization steps were reached, whichever occurred first.
