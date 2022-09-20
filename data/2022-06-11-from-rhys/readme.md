# File Index

Files in this folder were [sent by Rhys on Slack](https://ml-physics.slack.com/archives/DD8GBBRLN/p1654973643390109) on Jun 11, 2022.

They are needed to generate the classification histograms and moving MAE plots in the [Wren paper](https://arxiv.org/abs/2106.11132) and go together with the code in `mb_discovery/plots.py`.

## Format WBM material IDs

To convert `step_1_0`, `step_1_1` into `wbm-step-1-1`, `wbm-step-1-2`, etc., see `increment_wbm_material_id` in `process_wbm_cleaned.py`.

## WBM formation energies

All WBM `ComputedStructureEntries` have no energy corrections applied.

```py
import pandas as pd
from pymatgen.entries.computed_entries import ComputedStructureEntry

from mb_discovery import ROOT


df_wbm = pd.read_json(
    f"{ROOT}/data/2022-06-26-wbm-cses-and-initial-structures.json.gz"
).set_index("material_id")

CSEs = df_wbm.cse.map(ComputedStructureEntry.from_dict)

sum(cse.uncorrected_energy != cse.energy for cse in CSEs)
>>> 0
```
