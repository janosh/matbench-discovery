import pickle
from datetime import datetime

import pandas as pd


__author__ = "Janosh Riebesell"
__date__ = "2022-08-05"

today = f"{datetime.now():%Y-%m-%d}"

with open("block_0.p", "rb") as block_0:
    data = pickle.load(block_0)

with open("block_1.p", "rb") as block_1:
    data.update(pickle.load(block_1))


df = pd.DataFrame(data).T
df.index.name = "material_id"

# see https://github.com/materialsvirtuallab/m3gnet/issues/20#issuecomment-1206865218
# for why this is necessary
for structs in df.structure:
    for struct in structs:
        struct.lattice._pbc = (True, True, True)

df["final_structure"] = df.structure.map(lambda x: x[-1])
df["final_energy"] = df.energy.map(lambda x: x[-1])

df[["final_structure", "final_energy"]].reset_index().to_json(
    f"{today}-mp-2021-struct-energy.json.gz", default_handler=lambda x: x.as_dict()
)
