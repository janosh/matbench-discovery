"""Generate per-model hull-distance error box statistics for the site."""

from matbench_discovery import figs
from matbench_discovery.cli import complete_models
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.preds.discovery import df_each_err, df_preds

__author__ = "Janosh Riebesell"
__date__ = "2023-05-25"


test_subset = globals().get("test_subset", TestSubset.uniq_protos)

if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(MbdKey.uniq_proto)
    df_each_err = df_each_err.loc[df_preds.index]


# %%
box_models: list[dict[str, object]] = [
    {
        "key": model.key,
        "label": model.label,
        "quantiles": figs.round_list(
            df_each_err[model.label].quantile((0.05, 0.25, 0.5, 0.75, 0.95))
        ),
    }
    for model in complete_models()
]


# %%
figs.write_site_payload("box-hull-dist-errors", {"models": box_models})
