# %%
from diel_frontier.utils import decode_df_from_json

from ml_stability import ROOT


# %%
df = decode_df_from_json(f"{ROOT}/data/mp.tar.bz2", orient="split")


# %%
df = decode_df_from_json(f"{ROOT}/data/wbm_cleaned.json.gz")
