import os
from os.path import dirname


PKG_DIR = dirname(__file__)
ROOT = dirname(PKG_DIR)

os.makedirs(f"{PKG_DIR}/plots", exist_ok=True)
