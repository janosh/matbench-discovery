# %%
import pickle
from datetime import datetime

import pandas as pd
from aviary.wren.utils import get_aflow_label_from_spglib
from tqdm import tqdm

from matbench_discovery import ROOT

__author__ = "Janosh Riebesell"
__date__ = "2022-08-25"


today = f"{datetime.now():%Y-%m-%d}"


# %% block_{0,1}.p files downloaded from
# https://figshare.com/articles/dataset/MPF_2021_2_8/19470599
with open("block_0.p", "rb") as block_0:
    data = pickle.load(block_0)

with open("block_1.p", "rb") as block_1:
    data.update(pickle.load(block_1))


df_source = pd.DataFrame(data).T
df_source.index.name = "material_id"

# structures were originally pickled but following this comment reuploaded as CIFs.
# https://github.com/materialsvirtuallab/m3gnet/issues/20#issuecomment-1206865218
# see there for why this is necessary
for structs in df_source.structure:
    for struct in structs:
        struct.lattice._pbc = (True, True, True)


df = pd.DataFrame(index=df_source.index)
df["mp_structure"] = df.structure.map(lambda x: x[-1])
df["mp_energy"] = df.energy.map(lambda x: x[-1])


df["n_sites"] = df.mp_structure.map(len)
df["mp_energy_per_atom"] = df.mp_energy / df.n_sites

df["wyckoff"] = [get_aflow_label_from_spglib(x) for x in tqdm(df.mp_structure)]


df.reset_index().to_json(
    f"{ROOT}/dataset/{today}-m3gnet-trainset-mp-2021-struct-energy.json.gz",
    default_handler=lambda x: x.as_dict(),
)


# %% read processed M3GNet training data from disk
df = pd.read_json(
    f"{ROOT}/dataset/2022-08-25-m3gnet-trainset-mp-2021-struct-energy.json.gz"
).set_index("material_id")
