# %%
import itertools
import os
import subprocess
from collections.abc import Sequence
from datetime import date
from glob import glob

import numpy as np
import pandas as pd
import yaml
from pymatviz.enums import Key, eV_per_atom
from pymatviz.io import df_to_html, df_to_pdf
from pymatviz.utils import si_fmt
from sklearn.dummy import DummyClassifier

from matbench_discovery import DATA_DIR, PDF_FIGS, ROOT, SCRIPTS, SITE_FIGS
from matbench_discovery.data import DataFiles, df_wbm
from matbench_discovery.enums import MbdKey, Open, Targets
from matbench_discovery.metrics import stable_metrics
from matbench_discovery.models import MODEL_METADATA, model_is_compliant

try:
    from IPython.display import display
except ImportError:
    display = print

__author__ = "Janosh Riebesell"
__date__ = "2022-11-28"


# %% Create discovery metrics dataframes from MODEL_METADATA
def get_metrics_df(nested_keys: Sequence[str]) -> pd.DataFrame:
    """Extract metrics from MODEL_METADATA into a DataFrame.
    Returns a DataFrame with models as rows and metrics as columns."""
    out_dict = {}
    for model_name, metadata in MODEL_METADATA.items():
        metrics = None
        try:
            combined_metrics: dict[str, float] = {}
            for nested_key in nested_keys:
                metrics = metadata.get("metrics", {})
                for sub_key in nested_key.split("."):
                    metrics = metrics.get(sub_key, {})
                if not isinstance(metrics, dict):
                    break
                combined_metrics |= metrics
            if combined_metrics:
                out_dict[model_name] = combined_metrics

        except Exception as exc:
            exc.add_note(f"{model_name=} with {metrics=}")
            raise

    # Return DataFrame with models as rows instead of columns
    return pd.DataFrame.from_dict(out_dict, orient="index")


# Create DataFrames with models as rows
df_metrics = get_metrics_df(["discovery.full_test_set"]).sort_values(
    by=Key.f1, ascending=False
)
df_metrics_10k = get_metrics_df(["discovery.most_stable_10k"]).sort_values(
    by=Key.f1, ascending=False
)
df_metrics_uniq_protos = get_metrics_df(
    ["discovery.unique_prototypes", "phonons"]
).sort_values(by=Key.f1, ascending=False)
df_metrics_uniq_protos = df_metrics_uniq_protos.drop(
    columns=[MbdKey.missing_preds, MbdKey.missing_percent]
)

for df, title in (
    (df_metrics, "Metrics for Full Test Set"),
    (df_metrics_10k, "Metrics for 10k Most Stable Predictions"),
    (df_metrics_uniq_protos, "Metrics for unique non-MP prototypes"),
):
    df.attrs["title"] = title


# %%
name_map = {
    "MEGNet RS2RE": "MEGNet",
    "M3GNet→MEGNet": "M3GNet",
    "CHGNet→MEGNet": "CHGNet",
}
date_added_col = "Date Added"
model_name_col = "Model"

non_compliant_models = [
    key for key, meta in MODEL_METADATA.items() if not model_is_compliant(meta)
]

with open(f"{DATA_DIR}/training-sets.yml") as file:
    TRAINING_SETS = yaml.safe_load(file)

# Add model metadata to df_metrics(_10k, _uniq_protos)
for model in df_metrics_uniq_protos.index:
    if model == "Dummy":
        continue
    model_name = name_map.get(model, model)
    model_metadata = MODEL_METADATA.get(model_name, {})
    try:
        model_key = model_metadata.get("model_key", model_name)

        date_added = model_metadata.get("date_added", "")
        # long format date for tooltip, e.g. Monday, 28 November 2022
        title = f"{date.fromisoformat(date_added):%A, %d %B %Y}"
        df_metrics_uniq_protos.loc[model, date_added_col] = (
            f"<span {title=}>{date_added}</span>"
        )

        # Update targets column with full label in tooltip and data attribute
        model_targets = Targets[model_metadata[Key.targets]]
        tar_label = model_targets.label.replace(
            "<sub>", "<sub style='font-size: 0.8em;'>"
        )
        df_metrics_uniq_protos.loc[model, Key.targets.label] = (
            f'<span title="{model_targets.description}" '
            f'data-targets="{model_metadata[Key.targets]}">{tar_label}</span>'
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
        df_metrics_uniq_protos.loc[model, model_name_col] = (
            f"<span {html_attr_str}>{model}</span>"
        )

        if training_sets := model_metadata.get("training_set"):
            if isinstance(training_sets, dict | str):
                training_sets = [training_sets]

            n_structs_total = n_materials_total = 0
            dataset_urls, dataset_tooltip_lines = {}, []

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
                    dataset_tooltip_lines += [
                        f"{title}: {si_fmt(n_materials)} materials "
                        f"({si_fmt(n_structs)} structures)"
                    ]
                else:
                    dataset_tooltip_lines += [
                        f"{title}: {si_fmt(n_materials)} materials"
                    ]

            dataset_links = []
            for key, href in dataset_urls.items():
                if href:
                    rel = "noopener noreferrer"
                    dataset_links += [f"<a {href=} target='_blank' {rel=}>{key}</a>"]
                else:
                    dataset_links += [key]
            dataset_str = "+".join(dataset_links)
            new_line = "&#013;"  # line break that works in title attribute
            dataset_tooltip = (
                f"{new_line}&bull; ".join(["", *dataset_tooltip_lines])
                if len(dataset_tooltip_lines) > 1
                else ""
            )

            title = (
                f"{n_materials_total:,} materials in training set{new_line}"
                f"{dataset_tooltip}"
            )
            train_size_str = (
                f"<span {title=} data-sort-value={n_materials_total}>"
                f"{si_fmt(n_materials_total, fmt='.0f')} ({dataset_str})</span>"
            )

            if n_materials_total != n_structs_total:
                title = (
                    f"{n_materials_total:,} materials in training set "
                    f"({n_structs_total:,} structures counting all DFT relaxation "
                    f"frames per material){dataset_tooltip}"
                )
                train_size_str = (
                    f"<span {title=} data-sort-value={n_materials_total}>"
                    f"{si_fmt(n_materials_total, fmt='.0f')}"
                    f" <small>({si_fmt(n_structs_total, fmt='.1f')})</small>"
                    f" ({dataset_str})</span>"
                )

            df_metrics_uniq_protos.loc[model, Key.train_set.label] = train_size_str
        elif model == "Dummy":
            continue
        else:
            raise ValueError(
                f"Unknown {training_sets=} for {model=}\nwith {model_metadata=}"
            )

        model_params = model_metadata.get(Key.model_params, 0)
        n_estimators = model_metadata.get(Key.n_estimators, -1)
        title = f"{n_estimators:,} models in ensemble"
        n_estimators_str = (
            f" <small {title=}>(N={n_estimators})</small>" if n_estimators > 1 else ""
        )

        title = f"{model_params:,} trainable model parameters"
        formatted_params = si_fmt(model_params)
        df_metrics_uniq_protos.loc[
            model, Key.model_params.label.replace("eter", "")
        ] = (
            f'<span {title=} data-sort-value="{model_params}">{formatted_params}'
            f"</span>{n_estimators_str}"
        )

        for key in (MbdKey.openness, Key.train_task, Key.test_task):
            default = {MbdKey.openness: Open.OSOD}.get(key, pd.NA)
            df_metrics_uniq_protos.loc[model, key.label] = model_metadata.get(
                key, default
            )
    except Exception as exc:
        exc.add_note(f"{model=} with {model_metadata=}")
        raise

# assign this col to all tables
for df in (df_metrics, df_metrics_10k, df_metrics_uniq_protos):
    df[model_name_col] = df_metrics_uniq_protos[model_name_col]


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

    df_out.loc[col] = dummy_metrics | {model_name_col: col}


# %%
with open(f"{SCRIPTS}/metrics-which-is-better.yml") as file:
    better = yaml.safe_load(file)

R2_col = "R<sup>2</sup>"
kappa_srme_col = "κ<sub>SRME</sub>"
higher_is_better = {*better["higher_is_better"]} - {"R2"} | {R2_col}
lower_is_better = {*better["lower_is_better"]} | {kappa_srme_col}

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
]
show_cols = [
    *f"F1,DAF,Prec,Acc,TPR,TNR,MAE,RMSE,{R2_col},{kappa_srme_col}".split(","),
    *meta_cols,
]

for (label, df_met), show_non_compliant in itertools.product(
    (
        ("", df_metrics),
        ("-first-10k", df_metrics_10k),
        ("-uniq-protos", df_metrics_uniq_protos),
    ),
    (True, False),
):
    # abbreviate long column names
    df_met = df_met.rename(
        columns={
            "R2": R2_col,
            "Precision": "Prec",
            "Accuracy": "Acc",
            "κ_SRME": kappa_srme_col,
        }
    )
    # only keep columns we want to show
    df_table = df_met.filter([model_name_col, *show_cols])

    if make_uip_megnet_comparison:
        df_table = df_table[
            df_table.index.str.contains("MEGNet|CHGNet|M3GNet")
        ]  # |Dummy
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
            dict.fromkeys(df_table.select_dtypes(float), "{:.3g}"),
            na_rep="",  # render NaNs as empty string
        )
        .background_gradient(
            cmap="viridis", subset=list(higher_is_better & {*df_table}), axis="index"
        )
        .background_gradient(  # reverse color map if lower=better
            cmap="viridis_r", subset=list(lower_is_better & {*df_table}), axis="index"
        )
    )

    def tooltip_label(col: str) -> str:
        # add up/down arrows to indicate which metrics are better when higher/lower
        arrow_suffix = dict.fromkeys(higher_is_better, " ↑") | dict.fromkeys(
            lower_is_better, " ↓"
        )

        clf_suffix = (
            "of classifying thermodynamic stability (i.e. above/below 0K convex hull)"
        )
        reg_suffix = f"of predicting the convex hull distance {eV_per_atom}"
        tooltips_titles = {
            "Acc": f"accuracy {clf_suffix}",
            "DAF": "discovery acceleration factor",
            "F1": "harmonic mean of precision and recall",
            "Prec": f"precision {clf_suffix}",
            "TNR": f"true negative rate {clf_suffix}",
            "TPR": f"true positive rate {clf_suffix}",
            R2_col: "coefficient of determination",
            "MAE": f"mean absolute error {reg_suffix}",
            "RMSE": f"root mean squared error {reg_suffix}",
            "κ<sub>SRME</sub>": "symmetric relative mean error in predicted phonon "
            "mode contributions to thermal conductivity κ",
        }

        label = f"{col}{arrow_suffix.get(col, '')}"
        if col in tooltips_titles:
            title = tooltips_titles[col]
            return f"<span {title=}>{label}</span>"
        return label

    styler.relabel_index(
        [tooltip_label(col) for col in df_table], axis="columns"
    ).set_uuid("")

    # export model metrics as styled HTML table and Svelte component
    mae_col_idx = styler.columns.get_loc("MAE") + 1
    if kappa_srme_col in df_table:
        kappa_srme_col_idx = styler.columns.get_loc(kappa_srme_col) + 1
    else:
        kappa_srme_col_idx = 0

    # https://stackoverflow.com/a/38994837
    hide_scroll_bar = """
    table {
        scrollbar-width: none;  /* Firefox */
    }
    table::-webkit-scrollbar {
        display: none;  /* Safari and Chrome */
    }"""
    styler.set_sticky(axis="index")
    # Hide the original index since it's the same content same as model_name_col except
    # model_name_col also has HTML title attributes for hover tooltips
    styler.hide(axis="index")
    df_to_html(
        styler,
        file_path=f"{SITE_FIGS}/metrics-table{label}.svelte",
        inline_props="class='metrics'",
        # draw line between classification and regression metrics
        styles=f"#T_ :is(td, th):is(:nth-child({mae_col_idx}), "
        f":nth-child({kappa_srme_col_idx})) {{ border-left: 1px solid white; }}"
        f"{hide_scroll_bar}",
        sortable=True,
    )
    suffix = "" if show_non_compliant else "-only-compliant"
    non_compliant_idx = [*set(styler.index) & set(non_compliant_models)]
    try:
        os.environ["DYLD_FALLBACK_LIBRARY_PATH"] = "/opt/homebrew/lib"
        for pdf_path in (PDF_FIGS, f"{ROOT}/site/static/figs"):
            df_to_pdf(
                styler.hide([] if show_non_compliant else non_compliant_idx),
                f"{pdf_path}/metrics-table{label}{suffix}.pdf",
            )
    except (ImportError, RuntimeError, OSError) as exc:
        print(f"df_to_pdf failed: {exc}")

    display(styler.set_caption(df_met.attrs["title"]))


try:
    # convert PDFs in site/static/figs to SVGs (needed for "download as SVG" button)
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
