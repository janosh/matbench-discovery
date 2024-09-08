# %%
import itertools
import subprocess
from glob import glob

import numpy as np
import pandas as pd
import yaml
from pymatviz.enums import Key
from pymatviz.io import df_to_html_table, df_to_pdf
from pymatviz.utils import si_fmt
from sklearn.dummy import DummyClassifier

from matbench_discovery import DATA_DIR, PDF_FIGS, ROOT, SCRIPTS, SITE_FIGS
from matbench_discovery.data import DataFiles, df_wbm
from matbench_discovery.enums import MbdKey, Open
from matbench_discovery.metrics import stable_metrics
from matbench_discovery.models import MODEL_METADATA
from matbench_discovery.preds import df_metrics, df_metrics_10k, df_metrics_uniq_protos

try:
    from IPython.display import display
except ImportError:
    display = print

__author__ = "Janosh Riebesell"
__date__ = "2022-11-28"


# %%
name_map = {
    "MEGNet RS2RE": "MEGNet",
    "M3GNet→MEGNet": "M3GNet",
    "CHGNet→MEGNet": "CHGNet",
}
date_added_col = "Date Added"
df_metrics_uniq_protos.loc[Key.train_set.label] = df_metrics_uniq_protos.loc[
    date_added_col
] = ""

non_compliant_models = [
    key
    for key, meta in MODEL_METADATA.items()
    if meta.get("openness", Open.OSOD)
    != Open.OSOD  # TODO add `or uses_extra_training_data`
]

with open(f"{DATA_DIR}/training-sets.yml") as file:
    TRAINING_SETS = yaml.safe_load(file)

# add model metadata to df_metrics(_10k, _uniq_protos)
for model in df_metrics:
    model_name = name_map.get(model, model)
    model_metadata = MODEL_METADATA.get(model_name, {})
    model_key = model_metadata.get("model_key", model_name)

    df_metrics_uniq_protos.loc[date_added_col, model] = model_metadata.get(
        "date_added", ""
    )

    # Add model version as hover tooltip to model name
    model_version = model_metadata.get("model_version", "")
    compliant_css_cls = "non-compliant" if model in non_compliant_models else ""
    attrs = {
        "title": f"Version: {model_version}",
        "class": compliant_css_cls,
        "data-model-key": model_key,
    }
    html_attr_str = " ".join(f'{k}="{v}"' for k, v in attrs.items() if v)
    df_metrics_uniq_protos.loc[Key.model_name.label, model] = (
        f"<span {html_attr_str}>{model}</span>"
    )
    # assign this col to all tables
    for df in (df_metrics, df_metrics_10k, df_metrics_uniq_protos):
        df.loc[Key.model_name.label] = df_metrics_uniq_protos.loc[Key.model_name.label]

    if training_sets := model_metadata.get("training_set"):
        if isinstance(training_sets, dict | str):
            training_sets = [training_sets]

        n_structs_total = n_materials_total = 0
        dataset_urls, tooltip_lines = {}, []

        for train_set in training_sets:
            if isinstance(train_set, str) and train_set not in TRAINING_SETS:
                raise ValueError(f"Unknown training set {train_set=} for {model=}")
            key = train_set if isinstance(train_set, str) else ""
            dataset_info = TRAINING_SETS.get(key, train_set)
            n_structs = dataset_info["n_structures"]
            n_materials = dataset_info.get("n_materials", n_structs)

            n_structs_total += n_structs
            n_materials_total += n_materials

            title = dataset_info.get("title", key)
            dataset_urls[key or title] = dataset_info.get("url", "")

            if n_materials != n_structs:
                tooltip_lines += [
                    f"{title}: {si_fmt(n_materials)} materials ({si_fmt(n_structs)} "
                    "structures)"
                ]
            else:
                tooltip_lines += [f"{title}: {si_fmt(n_materials)} materials"]

        title = "Number of materials in training set"
        train_size_str = (
            f"<span {title=} data-sort-value={n_materials_total}>"
            f"{si_fmt(n_materials_total, fmt='.0f')}</span>"
        )

        if n_materials_total != n_structs_total:
            title = (
                "Number of materials in training set&#013;In parenthesis=number of "
                "structures, usually from multiple DFT relaxation frames per material"
            )
            train_size_str = (
                f"<span {title=} data-sort-value={n_materials_total}>"
                f"{si_fmt(n_materials_total, fmt='.0f')}</span>"
                f" <small {title=}> ({si_fmt(n_structs_total, fmt='.1f')})</small>"
            )

        dataset_links = []
        for key, href in dataset_urls.items():
            if href:
                dataset_links += [
                    f"<a {href=} target='_blank' rel='noopener noreferrer'>{key}</a>"
                ]
            else:
                dataset_links += [key]
        dataset_str = "+".join(dataset_links)

        tooltip = "&#013;".join(tooltip_lines)
        train_size_str = (
            f"<span title='{tooltip}'>{train_size_str} ({dataset_str})</span>"
        )

        df_metrics_uniq_protos.loc[Key.train_set.label, model] = train_size_str
    elif model == "Dummy":
        continue
    else:
        raise ValueError(
            f"Unknown {training_sets=} for {model=}\nwith {model_metadata=}"
        )

    model_params = model_metadata.get(Key.model_params, 0)
    n_estimators = model_metadata.get(Key.n_estimators, -1)
    title = "Number of models in ensemble"
    n_estimators_str = (
        f" <small {title=}>(N={n_estimators})</small>" if n_estimators > 1 else ""
    )

    title = "Number of trainable model parameters"
    formatted_params = si_fmt(model_params)
    df_metrics_uniq_protos.loc[Key.model_params.label.replace("eter", ""), model] = (
        f'<span {title=} data-sort-value="{model_params}">{formatted_params}'
        f"</span>{n_estimators_str}"
    )

    for key in (
        MbdKey.openness,
        Key.train_task,
        Key.test_task,
        Key.targets,
    ):
        default = {MbdKey.openness: Open.OSOD}.get(key, pd.NA)
        df_metrics_uniq_protos.loc[key.label, model] = model_metadata.get(key, default)


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

    df_out[col] = dummy_metrics | {Key.model_name.label: col}


# %%
with open(f"{SCRIPTS}/metrics-which-is-better.yml") as file:
    better = yaml.safe_load(file)

R2_col = "R<sup>2</sup>"
higher_is_better = {*better["higher_is_better"]} - {"R2"} | {R2_col}
lower_is_better = {*better["lower_is_better"]}

# if True, make metrics-table-megnet-uip-combos.(svelte|pdf) for SI
# if False, make metrics-table.(svelte|pdf) for main text
# when setting to True, uncomment the lines chgnet_megnet, m3gnet_megnet, megnet_rs2re
# in the Model enum!
make_uip_megnet_comparison = False

meta_cols = [
    Key.train_set.label,
    Key.model_params.label.replace("eter", ""),
    Key.model_type.label,
    Key.targets.label,
    date_added_col,
    Key.model_name.label,
]
show_cols = [*f"F1,DAF,Prec,Acc,TPR,TNR,MAE,RMSE,{R2_col}".split(","), *meta_cols]

for (label, df_met), show_non_compliant in itertools.product(
    (
        ("", df_metrics),
        ("-first-10k", df_metrics_10k),
        ("-uniq-protos", df_metrics_uniq_protos),
    ),
    (True, False),
):
    # abbreviate long column names
    df_met = df_met.rename(index={"R2": R2_col, "Precision": "Prec", "Accuracy": "Acc"})
    df_met.index.name = "Model"
    # only keep columns we want to show
    df_table = df_met.T.filter(show_cols)
    df_table = df_table.set_index(Key.model_name.label)
    df_table.index.name = None

    if make_uip_megnet_comparison:
        df_table = df_table.filter(regex="MEGNet|CHGNet|M3GNet")  # |Dummy
        label += "-uip-megnet-combos"
        if "M3GNet→MEGNet" not in df_table:
            print(
                "hint: for make_uip_megnet_comparison, uncomment the lines "
                "chgnet_megnet and m3gnet_megnet in the Model enum"
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
            cmap="viridis", subset=list(higher_is_better & {*df_table}), axis="index"
        )
        .background_gradient(  # reverse color map if lower=better
            cmap="viridis_r", subset=list(lower_is_better & {*df_table}), axis="index"
        )
    )
    # add up/down arrows to indicate which metrics are better when higher/lower
    arrow_suffix = dict.fromkeys(higher_is_better, " ↑") | dict.fromkeys(
        lower_is_better, " ↓"
    )
    styler.relabel_index(
        [f"{col}{arrow_suffix.get(col, '')}" for col in df_table], axis="columns"
    ).set_uuid("")

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
    suffix = "" if show_non_compliant else "-only-compliant"
    non_compliant_idx = [  # get HTML strings of non-compliant models in styler.index
        idx
        for idx in styler.index
        if any(f">{model_name}<" in idx for model_name in non_compliant_models)
    ]
    try:
        for pdf_path in (PDF_FIGS, f"{ROOT}/site/static/figs"):
            df_to_pdf(
                styler.hide([] if show_non_compliant else non_compliant_idx),
                f"{pdf_path}/metrics-table{label}{suffix}.pdf",
            )
    except (ImportError, RuntimeError) as exc:
        print(f"df_to_pdf failed: {exc}")

    display(styler.set_caption(df_met.attrs.get("title")))


try:
    # convert PDFs in site/static/figs to SVGs
    for pdf_path in glob(f"{ROOT}/site/static/figs/metrics-table*.pdf"):
        subprocess.run(
            ["pdf2svg", pdf_path, pdf_path.replace(".pdf", ".svg")], check=False
        )

    # svgo compress SVGs
    subprocess.run(["svgo", "--multipass", f"{ROOT}/site/static/figs"], check=False)
except FileNotFoundError:  # skip in CI where pdf2svg and svgo not installed
    pass


# %% PNG metrics table unused
if False:
    try:
        import dataframe_image

        png_metrics = f"{PDF_FIGS}/metrics-table.png"
        dataframe_image.export(styler, png_metrics, dpi=300)
    except ImportError:
        print("dataframe_image not installed, skipping png export")
