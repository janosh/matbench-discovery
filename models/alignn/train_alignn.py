"""Train a ALIGNN on target_col of data_path."""


# %%
from alignn.train_folder import train_for_folder

__author__ = "Philipp Benner"
__date__ = "2023-06-02"


# %%
train_for_folder(
    root_dir="data_train",
    config_name="alignn-config.json",
    keep_data_order=False,
    output_dir="data-train-result",
    epochs=1000,
    file_format="poscar",
)
