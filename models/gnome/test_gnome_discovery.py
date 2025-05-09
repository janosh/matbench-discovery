"""Only a placeholder script to inspect the model predictions. No way to reproduce
results since DeepMind did not release the model weights."""

# %%
import pymatviz as pmv

from matbench_discovery.enums import MbdKey, Model
from matbench_discovery.preds.discovery import df_preds

__author__ = "Janosh Riebesell"
__date__ = "2024-02-03"


# %%
df_preds[Model.gnome.label].hist(bins=100, figsize=(10, 10))


# %%
pmv.density_scatter(df=df_preds, x=MbdKey.e_form_dft, y=Model.gnome.label)
