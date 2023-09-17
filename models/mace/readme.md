## MACE formation energy predictions on WBM test set

This submission uses the [`2023-08-14-mace-yuan-trained-mptrj-04.model`](https://figshare.com/ndownloader/files/42374049) checkpoint trained by Yuan Chiang on the [MPtrj dataset](https://figshare.com/articles/dataset/23713842).
We initially tested the `2023-07-14-mace-universal-2-big-128-6.model` checkpoint trained on the much smaller [original M3GNet training set](https://figshare.com/articles/dataset/MPF_2021_2_8/19470599) which we received directly from Ilyes Batatia. MPtrj-trained MACE performed better and was used for the Matbench Discovery v1 submission.

MACE relaxed each test set structure until the maximum force in the training set dropped below 0.05 eV/Å or 500 optimization steps were reached, whichever occurred first.
