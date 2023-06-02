# %%

from alignn.train_folder import train_for_folder

"""
Train a ALIGNN on target_col of data_path.
"""

# %%

train_for_folder(
    root_dir='data_train',
    config_name='alignn_config.json',
    keep_data_order=False,
    output_dir='data_train_result',
    epochs=1000,
    file_format='poscar',
)
