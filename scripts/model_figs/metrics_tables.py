"""Compile metrics and total run times for all models and export them to JSON, a
pandas-styled HTML table and a plotly figure.
"""

# %%
import numpy as np
import pandas as pd
import yaml
from pymatviz.enums import Key
from pymatviz.io import df_to_html_table, df_to_pdf
from pymatviz.utils import si_fmt
from sklearn.dummy import DummyClassifier

from matbench_discovery import PDF_FIGS, SCRIPTS, SITE_FIGS
from matbench_discovery.data import DataFiles, df_wbm
from matbench_discovery.enums import MbdKey, ModelType, Open
from matbench_discovery.metrics import stable_metrics
from matbench_discovery.models import MODEL_METADATA
from matbench_discovery.preds import df_metrics, df_metrics_10k, df_metrics_uniq_protos

try:
    from IPython.display import display
except ImportError:
    display = print

__author__ = "Janosh Riebesell"
__date__ = "2022-11-28"


proprietaries = [
    key
    for key, meta in MODEL_METADATA.items()
    if meta.get("openness", Open.OSOD) != Open.OSOD
]


# %%
name_map = {
    "MEGNet RS2RE": "MEGNet",
    "M3GNet→MEGNet": "M3GNet",
    "CHGNet→MEGNet": "CHGNet",
}
df_met = df_metrics_uniq_protos
df_met.loc[Key.train_set.label] = ""

for model in df_metrics:
    model_name = name_map.get(model, model)
    if not (model_data := MODEL_METADATA.get(model_name)):
        continue
    n_structs = model_data["training_set"]["n_structures"]
    n_materials = model_data["training_set"].get("n_materials", n_structs)
    title = "Number of materials in training set"
    train_size_str = f"<span {title=}>{si_fmt(n_materials, fmt='.0f')}</span>"

    if n_materials != n_structs:
        title = "Number of materials (and structures) in training set"
        train_size_str = (
            f"<span {title=}>{si_fmt(n_materials, fmt='.0f')}</span> "
            f"<small {title=}>({si_fmt(n_structs, fmt='.0f')})</small>"
        )

    if train_url := model_data.get("training_set", {}).get("url"):
        train_size_str = (
            f"<a href='{train_url}' target='_blank' rel='noopener "
            f"noreferrer'>{train_size_str}</a>"
        )

    df_met.loc[Key.train_set.label, model] = train_size_str
    model_params = model_data.get(Key.model_params, "")
    n_estimators = model_data.get(Key.n_estimators, "")
    title = "Number of models in ensemble"
    n_estimators_str = (
        f" <small {title=}>(N={n_estimators})</small>" if n_estimators > 1 else ""
    )

    model_type = model_data.get(Key.model_type, "")
    title = ModelType.val_label_dict().get(model_type, "")
    df_met.loc[Key.model_type.label, model] = f"<span {title=}>{model_type}</span>"

    title = "Number of trainable model parameters"
    df_met.loc[Key.model_params.label.replace("eter", ""), model] = (
        f"<span {title=}>{si_fmt(model_params)}</span>{n_estimators_str}"
        if isinstance(model_params, int)
        else model_params
    )
    for key in (
        MbdKey.openness,
        Key.train_task,
        Key.test_task,
        Key.targets,
    ):
        default = {MbdKey.openness: Open.OSOD}.get(key, pd.NA)
        df_met.loc[key.label, model] = model_data.get(key, default)


# %% add dummy classifier results to df_metrics(_10k, _uniq_protos)
df_mp = pd.read_csv(DataFiles.mp_energies.path, index_col=0)

for df_in, df_out, col in (
    (df_wbm, df_metrics, "Dummy"),
    # "Dummy" for df_metrics_10k is still for the full test set, not dummy metrics on
    # only first 10k most stable predictions
    (df_wbm, df_metrics_10k, "Dummy"),
    (df_wbm.query(Key.uniq_proto), df_metrics_uniq_protos, "Dummy"),
):
    dummy_clf = DummyClassifier(strategy="stratified", random_state=0)
    dummy_clf.fit(np.zeros_like(df_mp[Key.each]), df_mp[Key.each] == 0)
    dummy_clf_preds = dummy_clf.predict(np.zeros(len(df_in)))

    each_true = df_in[MbdKey.each_true]
    dummy_metrics = stable_metrics(
        each_true, np.array([1, -1])[dummy_clf_preds.astype(int)], fillna=True
    )

    # important: regression metrics from dummy_clf are meaningless, we overwrite them
    # with correct values here. don't remove!
    dummy_metrics[Key.daf] = 1
    dummy_metrics["R2"] = 0
    dummy_metrics["MAE"] = (each_true - each_true.mean()).abs().mean()
    dummy_metrics["RMSE"] = ((each_true - each_true.mean()) ** 2).mean() ** 0.5

    df_out[col] = dummy_metrics


# %%
with open(f"{SCRIPTS}/metrics-which-is-better.yml") as file:
    better = yaml.safe_load(file)

R2_col = "R<sup>2</sup>"
higher_is_better = {*better["higher_is_better"]} - {"R2"} | {R2_col}
lower_is_better = {*better["lower_is_better"]}

# if True, make metrics-table-megnet-uip-combos.(svelte|pdf) for SI
# if False, make metrics-table.(svelte|pdf) for main text
# when setting to True, uncomment the lines chgnet_megnet, m3gnet_megnet, megnet_rs2re
# in PredFiles!
make_uip_megnet_comparison = False
meta_cols = [
    Key.train_set.label,
    Key.model_params.label.replace("eter", ""),
    Key.model_type.label,
    Key.targets.label,
]
show_cols = [*f"F1,DAF,Prec,Acc,TPR,TNR,MAE,RMSE,{R2_col}".split(","), *meta_cols]

for label, df_met in (
    ("", df_metrics),
    ("-uniq-protos", df_metrics_uniq_protos),
    ("-first-10k", df_metrics_10k),
):
    # abbreviate long column names
    df_met = df_met.rename(index={"R2": R2_col, "Precision": "Prec", "Accuracy": "Acc"})
    df_met.index.name = "Model"
    # only keep columns we want to show
    df_table = df_met.T.filter(show_cols)

    if make_uip_megnet_comparison:
        df_table = df_table.filter(regex="MEGNet|CHGNet|M3GNet")  # |Dummy
        label += "-uip-megnet-combos"
        if "M3GNet→MEGNet" not in df_table:
            print(
                "hint: for make_uip_megnet_comparison, uncomment the lines "
                "chgnet_megnet and m3gnet_megnet in PredFiles"
            )

    if "-first-10k" in label:
        # hide redundant metrics for first 10k preds (all TPR = 1, TNR = 0)
        df_table = df_table.drop(["TPR", "TNR"], axis="columns")
    if label != "-uniq-protos":  # only show training size and model type once
        df_table = df_table.drop(meta_cols, axis="columns", errors="ignore")

    styler = (
        df_table.style.format(
            # render integers without decimal places
            dict.fromkeys("TP FN FP TN".split(), "{:,.0f}"),
            precision=3,  # render floats with 2 decimals
            na_rep="",  # render NaNs as empty string
        )
        .background_gradient(
            cmap="viridis", subset=list(higher_is_better & {*df_table}), axis=0
        )
        .background_gradient(  # reverse color map if lower=better
            cmap="viridis_r", subset=list(lower_is_better & {*df_table}), axis=0
        )
    )
    # add up/down arrows to indicate which metrics are better when higher/lower
    arrow_suffix = dict.fromkeys(higher_is_better, " ↑") | dict.fromkeys(
        lower_is_better, " ↓"
    )
    styler.relabel_index(
        [f"{col}{arrow_suffix.get(col, '')}" for col in df_table], axis="columns"
    ).set_uuid("")

    # add CSS class 'proprietary' to cells of proprietary models (openness != OSOD)
    styler.set_td_classes(
        df_table.T.assign(**dict.fromkeys(proprietaries, "proprietary"))[
            proprietaries
        ].T
    )

    # export model metrics as styled HTML table and Svelte component
    # get index of MAE column
    mae_col_idx = styler.columns.get_loc("MAE")
    col_selector = f"#T_ :is(td, th):nth-child({mae_col_idx + 2})"
    # https://stackoverflow.com/a/38994837
    hide_scroll_bar = """
    table {
        scrollbar-width: none;  /* Firefox */
    }
    table::-webkit-scrollbar {
        display: none;  /* Safari and Chrome */
    }"""
    df_to_html_table(
        styler,
        file_path=f"{SITE_FIGS}/metrics-table{label}.svelte",
        inline_props="class='metrics'",
        # draw line between classification and regression metrics
        styles=f"{col_selector} {{ border-left: 1px solid white; }}{hide_scroll_bar}",
    )
    try:
        df_to_pdf(styler, f"{PDF_FIGS}/metrics-table{label}.pdf")
    except (ImportError, RuntimeError) as exc:
        print(f"df_to_pdf failed: {exc}")

    display(styler.set_caption(df_met.attrs.get("title")))


# %% PNG metrics table unused
if False:
    try:
        import dataframe_image

        png_metrics = f"{PDF_FIGS}/metrics-table.png"
        dataframe_image.export(styler, png_metrics, dpi=300)
    except ImportError:
        print("dataframe_image not installed, skipping png export")
