"""This script can still be used to create the materials discovery metrics tables
locally but is no longer relevant for the online leaderboard which directly reads
metrics from the model YAML files and constructs the table on the fly.

It is still useful for creating the PDF and SVG versions of the table for download from
the site.
"""

# %%
import traceback
from datetime import date

import numpy as np
import pandas as pd
import pymatviz as pmv
import yaml
from pymatviz.enums import Key, eV_per_atom
from pymatviz.utils import si_fmt
from sklearn.dummy import DummyClassifier

from matbench_discovery import PKG_DIR, STABILITY_THRESHOLD
from matbench_discovery.cli import cli_args
from matbench_discovery.data import DATASETS, df_wbm
from matbench_discovery.enums import DataFiles, MbdKey, Model, Open, Targets, TestSubset
from matbench_discovery.metrics import discovery

__author__ = "Janosh Riebesell"
__date__ = "2022-11-28"


# %%
if __name__ == "__main__":
    for model in cli_args.models:
        if not model.is_complete:
            print(f"\nSkipping {model.label}: incomplete discovery metrics")
    cli_args.models = [model for model in cli_args.models if model.is_complete]
    if not cli_args.models:
        raise SystemExit(0)

    from matbench_discovery.preds.discovery import df_each_pred, df_preds

    uniq_protos_idx = df_wbm.query(MbdKey.uniq_proto).index
    # dummy discovery rate of stable crystals when selecting randomly from the unique
    # prototype subset, used to compute the discovery acceleration factor (DAF)
    uniq_proto_prevalence = (
        df_wbm.query(MbdKey.uniq_proto)[MbdKey.each_true] <= STABILITY_THRESHOLD
    ).mean()

    for model in cli_args.models:
        try:
            print(f"\nProcessing {model.label}...")
            model_preds = df_preds[model.label]
            each_true = df_preds[MbdKey.each_true]
            each_pred = df_each_pred[model.label]
            each_pred_uniq_protos = each_pred.loc[uniq_protos_idx]
            # 10k most stable = lowest *predicted hull distance* (each_pred) within
            # the unique prototype subset, not lowest raw formation energy
            most_stable_10k = each_pred_uniq_protos.nsmallest(10_000)

            full_metrics = discovery.stable_metrics(each_true, each_pred, fillna=True)
            uniq_proto_metrics = discovery.stable_metrics(
                each_true.loc[uniq_protos_idx], each_pred_uniq_protos, fillna=True
            )
            stable_10k_metrics = discovery.stable_metrics(
                each_true.loc[most_stable_10k.index], most_stable_10k, fillna=True
            )
            # DAF on these subsets is relative to the uniq-proto prevalence
            for metrics in (uniq_proto_metrics, stable_10k_metrics):
                metrics[str(Key.daf.symbol)] = (
                    metrics["Precision"] / uniq_proto_prevalence
                )

            for test_subset, (metrics, subset_idx) in {
                TestSubset.full_test_set: (full_metrics, slice(None)),
                TestSubset.uniq_protos: (uniq_proto_metrics, uniq_protos_idx),
                TestSubset.most_stable_10k: (
                    stable_10k_metrics,
                    most_stable_10k.index,
                ),
            }.items():
                # cast all values to float (incl. TP/FP/TN/FN counts) to match the
                # number format of previous YAML writes
                rounded_metrics: dict[str, str | float] = {
                    key: round(float(val), 3) for key, val in metrics.items()
                }
                discovery.write_metrics_to_yaml(
                    model, rounded_metrics, model_preds.loc[subset_idx], test_subset
                )
                print(f"\tUpdated discovery metrics for {test_subset}")
        except (ValueError, OSError, KeyError):
            print(f"\tError processing {model.label}: {traceback.format_exc()}")
            continue

    if not pmv.IS_IPYTHON:
        raise SystemExit(0)


# %%
date_added_col = "Date Added"
model_name_col = "Model"

# Add model metadata to df_metrics(_10k|_uniq_protos)
models = discovery.df_metrics_uniq_protos.columns
for model_label in models:
    if model_label == "Dummy":
        continue
    # Find the model enum entry by label
    model = Model.from_label(model_label)
    if not model.is_complete:
        continue

    try:
        model_key = model.metadata["model_key"]

        date_added = model.metadata.get("date_added", "")
        # long format date for tooltip, e.g. Monday, 28 November 2022
        title = f"{date.fromisoformat(date_added):%A, %d %B %Y}"
        discovery.df_metrics_uniq_protos.loc[date_added_col, model] = (
            f"<span {title=}>{date_added}</span>"
        )

        # Update targets column with full label in tooltip and data attribute
        model_targets = Targets[model.metadata[Key.targets]]
        tar_label = model_targets.label.replace(
            "<sub>", "<sub style='font-size: 0.8em;'>"
        )
        discovery.df_metrics_uniq_protos.loc[Key.targets.label, model] = (
            f'<span title="{model_targets.description}" '
            f'data-targets="{model.metadata[Key.targets]}">{tar_label}</span>'
        )

        # Add model version as hover tooltip to model name
        model_version = model.metadata.get("model_version", "")
        attrs = {
            "title": f"Version: {model_version}",
            "data-model-key": model_key,
        }
        html_attr_str = " ".join(f'{k}="{v}"' for k, v in attrs.items() if v)
        discovery.df_metrics_uniq_protos.loc[model_name_col, model] = (
            f"<span {html_attr_str}>{model}</span>"
        )

        if training_sets := model.metadata.get("training_set"):
            if isinstance(training_sets, dict | str):
                training_sets = [training_sets]

            n_structs_total = n_materials_total = 0
            dataset_urls, dataset_tooltip_lines = {}, []

            for train_set in training_sets:
                if isinstance(train_set, str) and train_set not in DATASETS:
                    raise ValueError(f"Unknown training set {train_set=} for {model=}")
                key = train_set if isinstance(train_set, str) else ""
                dataset_info = DATASETS.get(key, train_set)
                n_structs = dataset_info["n_structures"]
                n_materials = dataset_info.get("n_materials", n_structs)

                n_structs_total += n_structs
                n_materials_total += n_materials

                title = dataset_info.get("title", key)
                dataset_urls[key or title] = dataset_info.get("download_url", "")

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

            title = f"{n_materials_total:,} materials in training set{new_line}{dataset_tooltip}"  # noqa: E501
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

            discovery.df_metrics_uniq_protos.loc[Key.train_set.label, model] = (
                train_size_str
            )
        elif model == "Dummy":
            continue
        else:
            raise ValueError(
                f"Unknown {training_sets=} for {model=}\nwith {model.metadata=}"
            )

        model_params = model.metadata.get(Key.model_params, 0)
        n_estimators = model.metadata.get(Key.n_estimators, -1)
        title = f"{n_estimators:,} models in ensemble"
        n_estimators_str = (
            f" <small {title=}>(N={n_estimators})</small>" if n_estimators > 1 else ""
        )

        title = f"{model_params:,} trainable model parameters"
        formatted_params = si_fmt(model_params)
        discovery.df_metrics_uniq_protos.loc[
            Key.model_params.label.replace("eter", ""), model
        ] = (
            f'<span {title=} data-sort-value="{model_params}">{formatted_params}'
            f"</span>{n_estimators_str}"
        )

        for key in (MbdKey.openness, Key.train_task, Key.test_task):
            default = Open.OSOD if key == MbdKey.openness else pd.NA
            discovery.df_metrics_uniq_protos.loc[key.label, model] = model.metadata.get(
                key, default
            )
    except Exception as exc:
        exc.add_note(f"{model=} with {model.metadata=}")
        raise

# assign this col to all tables
for df in (
    discovery.df_metrics,
    discovery.df_metrics_10k,
    discovery.df_metrics_uniq_protos,
):
    df.loc[model_name_col, :] = discovery.df_metrics_uniq_protos.loc[model_name_col, :]


# %% add dummy classifier results to discovery.df_metrics(_10k, _uniq_protos)
df_mp = pd.read_csv(DataFiles.mp_energies.path, index_col=0)

for df_in, df_out, col in (
    (df_wbm, discovery.df_metrics, "Dummy"),
    # "Dummy" for df_metrics_10k is still for the full test set, not dummy metrics on
    # only first 10k most stable predictions
    (df_wbm, discovery.df_metrics_10k, "Dummy"),
    (df_wbm.query(MbdKey.uniq_proto), discovery.df_metrics_uniq_protos, "Dummy"),
):
    dummy_clf = DummyClassifier(strategy="stratified", random_state=0)
    dummy_clf.fit(
        np.zeros_like(df_mp["energy_above_hull"]), df_mp["energy_above_hull"] == 0
    )
    dummy_clf_preds = dummy_clf.predict(np.zeros(len(df_in)))

    each_true = df_in[MbdKey.each_true]
    dummy_metrics = discovery.stable_metrics(
        each_true, np.array([1, -1])[dummy_clf_preds.astype(int)], fillna=True
    )

    # important: regression metrics from dummy_clf are meaningless, we overwrite them
    # with correct values here. don't remove!
    dummy_metrics[str(Key.daf.symbol)] = 1
    dummy_metrics["R2"] = 0
    dummy_metrics["MAE"] = (each_true - each_true.mean()).abs().mean()
    dummy_metrics["RMSE"] = ((each_true - each_true.mean()) ** 2).mean() ** 0.5

    df_out.loc[:, col] = dummy_metrics | {model_name_col: col}


# %%
with open(f"{PKG_DIR}/modeling-tasks.yml", encoding="utf-8") as file:
    discovery_metrics = yaml.safe_load(file)["discovery"]["metrics"]

R2_col = "R<sup>2</sup>"
kappa_srme_col = "κ<sub>SRME</sub>"
higher_is_better = {*discovery_metrics["higher_is_better"]} - {"R2"} | {R2_col}
lower_is_better = {*discovery_metrics["lower_is_better"]} | {kappa_srme_col}

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
    *f"F1,DAF,Prec,Acc,TPR,TNR,MAE,RMSE,{R2_col}".split(","),
    kappa_srme_col,
    *meta_cols,
]

for label, df_met in (
    ("", discovery.df_metrics),
    ("-first-10k", discovery.df_metrics_10k),
    ("-uniq-protos", discovery.df_metrics_uniq_protos),
):
    out_label = label
    df_show = df_met.copy().T.rename_axis(index=model_name_col)
    if "kappa_103" in df_show:
        df_kappa = pd.json_normalize(df_show["kappa_103"])
        df_kappa.index = df_show.index
        df_show["κ_SRME"] = df_kappa["κ_SRME"]

    # abbreviate long column names
    df_show = df_show.rename(
        columns={
            "R2": R2_col,
            "Precision": "Prec",
            "Accuracy": "Acc",
            "κ_SRME": kappa_srme_col,
        }
    )
    # only keep columns we want to show
    df_table = df_show.filter([model_name_col, *show_cols])
    df_table = df_table.convert_dtypes()
    if make_uip_megnet_comparison:
        df_table = df_table[
            df_table.index.str.contains("MEGNet|CHGNet|M3GNet")
        ]  # |Dummy
        out_label += "-uip-megnet-combos"
        if "M3GNet→MEGNet" not in df_table:
            print(
                "hint: for make_uip_megnet_comparison, uncomment the lines "
                "chgnet_megnet and m3gnet_megnet in the Model enum"
            )

    if "-first-10k" in out_label:
        # hide redundant metrics for first 10k preds (all TPR = 1, TNR = 0)
        df_table = df_table.drop(["TPR", "TNR"], axis="columns")
    if out_label != "-uniq-protos":  # only show training size and model type once
        df_table = df_table.drop(meta_cols, axis="columns", errors="ignore")

    styler = df_table.style.format(
        dict.fromkeys(df_table.select_dtypes(float), "{:,.3f}"),  # use for manuscript
        na_rep="",  # render NaNs as empty string
    )
    styler = styler.background_gradient(  # ty: ignore[unresolved-attribute]  # pandas types .format() as StylerRenderer
        cmap="viridis", subset=list(higher_is_better & {*df_table}), axis="index"
    ).background_gradient(  # reverse color map if lower=better
        cmap="viridis_r", subset=list(lower_is_better & {*df_table}), axis="index"
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

    styler.set_sticky(axis="index")
    # Hide the original index since it's the same content same as model_name_col except
    # model_name_col also has HTML title attributes for hover tooltips
    styler.hide(axis="index")
    if pmv.IS_IPYTHON:
        from IPython.display import display

        display(styler.set_caption(df_show.attrs["title"]))
