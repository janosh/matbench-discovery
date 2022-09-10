# %%
from datetime import datetime
from glob import glob

import pandas as pd

__author__ = "Janosh Riebesell"
__date__ = "2022-08-03"

"""
This module loads all the data files in this directory as pandas dataframes and writes
the output of their describe() method to markdown file.
"""

today = f"{datetime.now():%Y-%m-%d}"


# %%
with open("data-describe.md", "w") as md_file:
    md_file.write("# Data File Stats\n\n")

    for filename in glob("*.csv"):
        df = pd.read_csv(filename)
        md_file.write(f"## `{filename}`\n\n")
        md_file.write(f"{df.shape = }\n\n")
        md_file.write(f"columns = {', '.join(df)}\n\n")
