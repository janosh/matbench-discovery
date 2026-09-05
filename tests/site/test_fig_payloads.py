"""Shape guard for the committed data-only figure payloads in site/src/figs.

The Svelte pages type these payloads via the ambient module declarations in
site/src/figs/payloads.d.ts, so shape drift between the Python exporters and the site
only surfaces as broken figures. These tests mirror payloads.d.ts and fail fast.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import os
from functools import partial
from typing import Any

import pandas as pd
import pytest

from matbench_discovery import ROOT, SITE_DIR, SITE_FIG_DATA, STABILITY_THRESHOLD, figs
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.enums import MbdKey, Model, TestSubset
from matbench_discovery.metrics.discovery import stable_metrics
from scripts.model_figs.kappa_103_analysis import row_flag


def reject_json_constant(const: str) -> None:
    """Reject NaN/Infinity literals in JSON figure payloads."""
    raise ValueError(f"non-finite JSON constant {const!r}")


def load_payload(name: str) -> dict[str, Any]:
    """Load a committed payload by stem: gzipped ``<name>.json.gz`` or ``<name>.jsonl``
    (reassembled like the site's jsonl Vite plugin). The exporter writes
    ``allow_nan=False``; the aggregate path also rejects NaN/Infinity to catch drift.
    """
    path = f"{SITE_FIG_DATA}/{name}.json.gz"
    if os.path.isfile(path):
        with gzip.open(path) as file:
            return json.load(file, parse_constant=reject_json_constant)
    return figs.read_jsonl_payload(f"{SITE_FIG_DATA}/{name}.jsonl")


def assert_num_list(values: object, *, length: int | None = None) -> None:
    """Assert a list of finite numbers (None allowed for gaps), optionally sized."""
    assert isinstance(values, list), f"expected list, got {type(values)=}"
    assert values, "expected non-empty list"
    # None marks a gap; every other entry must be a finite number. The exporter writes
    # allow_nan=False, so NaN/Infinity here means corruption - reject it on the .jsonl
    # path too (load_payload can't apply json's parse_constant guard line-by-line).
    assert all(
        val is None or (isinstance(val, (int, float)) and math.isfinite(val))
        for val in values
    )
    if length is not None:
        assert len(values) == length, f"expected {length=}, got {len(values)}"


def assert_xy(obj: dict[str, Any], *, bar_width: bool = False) -> None:
    """Assert an XY (or HistBins) payload fragment with matching lengths."""
    assert_num_list(obj["x"])
    assert_num_list(obj["y"], length=len(obj["x"]))
    if bar_width:
        assert isinstance(obj["bar_width"], (int, float))


def assert_models(
    payload: dict[str, Any], *keys: str, n_min: int = 2
) -> list[dict[str, Any]]:
    """Assert a {models: [...]} payload where each entry has (at least) ``keys``."""
    models = payload["models"]
    assert len(models) >= n_min, f"expected >= {n_min} models, got {len(models)}"
    for entry in models:
        for key in keys:
            assert key in entry, f"model {entry['label']!r} missing {key!r}"
    return models


def assert_model_keys(model: dict[str, Any], *derived_keys: str) -> None:
    """Assert a model record has exactly its identity and declared derived fields."""
    assert set(model) == {
        "model_key",
        "label",
        "input_artifacts",
        *derived_keys,
    }


def payload_model_keys(name: str) -> set[str]:
    """Set of immutable model keys in one payload."""
    return {model["model_key"] for model in load_payload(name)["models"]}


def test_no_orphan_payloads() -> None:
    """Every committed payload has a shape test below (and thus a consumer page)."""
    aggregates: set[str] = set()  # gzipped <name>.json.gz (static payloads)
    jsonl: set[str] = set()  # line-delimited <name>.jsonl (multi-model payloads)
    for entry in os.listdir(SITE_FIG_DATA):
        if entry.endswith(".json.gz"):
            aggregates.add(entry.removesuffix(".json.gz"))
        elif entry.endswith(".jsonl"):
            jsonl.add(entry.removesuffix(".jsonl"))
    # a payload is EITHER a gzipped aggregate or a .jsonl, never both (a stray leftover
    # aggregate beside a .jsonl would reintroduce the merge conflict .jsonl avoids)
    assert not (both := aggregates & jsonl), f"both aggregate and jsonl: {both}"
    on_disk = aggregates | jsonl
    assert on_disk == set(EXPECTED_PAYLOADS), (
        f"unexpected={on_disk - set(EXPECTED_PAYLOADS)}, "
        f"missing={set(EXPECTED_PAYLOADS) - on_disk}"
    )


def test_per_element_each_errors_payload() -> None:
    """Each per-element-each-errors column maps a key to finite per-element values."""
    path = f"{SITE_DIR}/routes/models/per-element-each-errors.jsonl"
    payload = figs.read_jsonl_payload(path)
    assert_num_list(list(payload["mp_occurrences"].values()))
    assert_num_list(list(payload["test_set_standard_deviation"].values()))
    columns = payload["models"]
    assert len(columns) > 10, f"expected >10 columns, got {len(columns)}"
    for column in columns:
        assert isinstance(column["model_key"], str)
        assert_num_list(list(column["values"].values()))


def check_box_hull_dist_errors() -> None:
    for model in assert_models(load_payload("box-hull-dist-errors")):
        assert_num_list(model["quantiles"], length=5)  # q05, q25, median, q75, q95


def check_cumulative_precision_recall() -> None:
    payload = load_payload("cumulative-precision-recall")
    assert payload["n_stable"] > 10_000
    for model in assert_models(payload):
        assert_num_list(model["x"])
        assert_num_list(model["precision"], length=len(model["x"]))
        assert_num_list(model["recall"], length=len(model["x"]))
        assert_num_list(model["end"], length=3)


def check_roc_models() -> None:
    for model in assert_models(load_payload("roc-models"), "auc"):
        assert 0.5 < model["auc"] <= 1
        assert_num_list(model["fpr"])
        assert_num_list(model["tpr"], length=len(model["fpr"]))


def check_rolling_mae() -> None:
    payload = load_payload("rolling-mae-vs-hull-dist")
    assert_num_list(payload["x"])
    for model in assert_models(payload):
        assert_num_list(model["y"], length=len(payload["x"]))
    assert_xy(payload["density"])


def check_hist_clf() -> None:
    payload = load_payload("hist-clf-pred-hull-dist")
    assert_num_list(payload["bin_centers"])
    n_bins = len(payload["bin_centers"])
    for model in assert_models(payload, "f1"):
        for clf in ("tp", "fn", "fp", "tn"):
            assert_num_list(model[clf], length=n_bins)

    record = next(
        (
            model
            for model in payload["models"]
            if model["f1"]
            != Model.from_ref(model["model_key"]).metrics["discovery"][
                TestSubset.uniq_protos.value
            ]["F1"]
        ),
        None,
    )
    # payload F1 values all match the registry fixtures: the shape checks above already
    # passed, so there is no mismatched record left to recompute against predictions
    if record is None:
        return
    test_model = Model.from_ref(record["model_key"])
    df_preds = load_df_wbm_with_preds(
        models=[test_model], subset=TestSubset.uniq_protos
    )
    each_pred = (
        df_preds[MbdKey.each_true]
        + df_preds[test_model.key]
        - df_preds[MbdKey.e_form_dft]
    )
    expected_f1 = round(
        stable_metrics(
            df_preds[MbdKey.each_true],
            each_pred,
            stability_threshold=STABILITY_THRESHOLD,
        )["F1"],
        4,
    )
    assert record["f1"] == expected_f1


def check_element_prevalence() -> None:
    payload = load_payload("element-prevalence-vs-error")
    elements = payload["elements"]
    assert all(isinstance(el, str) for el in elements)
    assert_num_list(payload["occurrences"], length=len(elements))
    for model in assert_models(payload):
        assert_model_keys(model, "y")
        assert_num_list(model["y"], length=len(elements))


def check_scatter_largest_fp_diff() -> None:
    payload = load_payload("scatter-largest-fp-diff-each-error")
    assert_num_list(payload["fp_diff"])
    for model in assert_models(payload, "mae"):
        assert_model_keys(model, "mae", "y")
        assert_num_list(model["y"], length=len(payload["fp_diff"]))


def check_hist_largest_each_errors() -> None:
    for model in assert_models(load_payload("hist-largest-each-errors-fp-diff")):
        assert_model_keys(model, "err_min", "err_max")
        assert_xy(model["err_min"], bar_width=True)
        assert_xy(model["err_max"], bar_width=True)


def check_hist_wbm_e_form_per_atom() -> None:
    assert_xy(load_payload("hist-wbm-e-form-per-atom"), bar_width=True)


def check_hist_wbm_hull_dist() -> None:
    payload = load_payload("hist-wbm-hull-dist")
    assert_xy(payload["stable"])
    assert_xy(payload["unstable"])
    assert isinstance(payload["bar_width"], float)
    assert payload["std"] > 0


def check_spacegroup_sunbursts() -> None:
    payload = load_payload("spacegroup-sunbursts")
    for key in ("mp", "wbm"):
        sunburst = payload[key]
        labels, parents = sunburst["labels"], sunburst["parents"]
        values, ids = sunburst["values"], sunburst["ids"]
        n_nodes = len(labels)
        assert len(parents) == len(values) == len(ids) == n_nodes > 0
        assert parents.count("") == 7  # the 7 crystal systems are the root nodes
        assert all(labels)
        assert all(val > 0 for val in values)
        assert len(set(ids)) == n_nodes  # ids unique


def check_arity_hist() -> None:
    payload = load_payload("mp-vs-mp-trj-vs-wbm-arity-hist")
    assert len(payload["datasets"]) == 3  # MP, MPtrj, WBM
    for dataset in payload["datasets"]:
        assert_xy(dataset)


def check_mp_trj_hists() -> None:
    payload = load_payload("mp-trj-hists")
    for key in ("e-form", "forces", "stresses", "magmoms", "n-sites"):
        assert_xy(payload[key], bar_width=True)
    n_sites = payload["n-sites"]
    assert_num_list(n_sites["cumulative"], length=len(n_sites["x"]))
    assert n_sites["cumulative"][-1] == pytest.approx(1, abs=1e-4)


def check_mp_elemental_ref_energies() -> None:
    assert_xy(load_payload("mp-elemental-ref-energies"))


def check_element_counts() -> None:
    payload = load_payload("element-counts-mp-vs-wbm")
    for variant in ("raw", "normalized"):
        assert len(payload[variant]) == 2  # WBM + MP
        for series in payload[variant]:
            assert all(isinstance(symbol, str) for symbol in series["x"])
            assert_num_list(series["y"], length=len(series["x"]))


def check_spg_sankeys() -> None:
    for model in assert_models(load_payload("spg-sankeys")):
        labels = model["labels"]
        source, target, value = model["source"], model["target"], model["value"]
        n_links = len(source)
        assert len(target) == len(value) == n_links > 0
        assert all(labels)
        for idx in range(n_links):
            assert 0 <= source[idx] < len(labels)
            assert 0 <= target[idx] < len(labels)
            assert value[idx] > 0


def check_kappa_103_analysis() -> None:
    payload = load_payload("kappa-103-analysis")
    material_ids = payload["material_ids"]
    n_materials = len(material_ids)
    assert n_materials == 103
    assert all(isinstance(mid, str) and mid for mid in material_ids)
    assert all(isinstance(formula, str) for formula in payload["formulas"])
    assert len(payload["formulas"]) == n_materials
    assert len(payload["spg_nums"]) == n_materials
    assert all(spg is None or 1 <= spg <= 230 for spg in payload["spg_nums"])
    assert_num_list(payload["kappa_dft"], length=n_materials)
    for model in assert_models(payload, "freq_w1_mean", "freq_pairs"):
        for field in ("kappa_ml", "srme", "freq_w1"):
            assert_num_list(model[field], length=n_materials)
        assert all(val is None or 0 <= val <= 2 for val in model["srme"])
        for field in ("imag_modes", "broken_sym", "max_steps"):
            flags = model[field]
            assert len(flags) == n_materials
            assert all(flag is None or isinstance(flag, bool) for flag in flags)
        # freq_pairs may be empty for models whose phonon runs all failed
        pairs = model["freq_pairs"]
        if pairs["dft"] or pairs["ml"]:
            assert_num_list(pairs["dft"])
            assert_num_list(pairs["ml"], length=len(pairs["dft"]))
        else:
            assert model["freq_w1_mean"] is None


@pytest.mark.parametrize(
    "value, expected",
    [
        (None, None),
        (float("nan"), None),
        (pd.NA, None),
        (1, True),
        (0, False),
    ],
)
def test_kappa_103_analysis_row_flag(value: object, expected: bool | None) -> None:
    """Failure flags treat missing values as absent and present values as bools."""
    assert row_flag(pd.Series({"flag": value}), "flag") is expected


def check_xy_models(name: str, *stat_keys: str) -> None:
    """Generic check: per-model x/y series plus data-derived stat fields."""
    for model in assert_models(load_payload(name), *stat_keys):
        if name == "scatter-largest-each-errors-fp-diff":
            assert_model_keys(model, "mae", "x", "y")
        assert_xy(model)


# payloads whose models are plain x/y series with the given stat fields
XY_MODEL_STATS = {
    "scatter-largest-each-errors-fp-diff": ("mae",),
    "struct-rmsd-cdf": ("auc",),
    "sym-ops-diff-bar": ("sigma",),
}

EXPECTED_PAYLOADS = {
    **{
        name: partial(check_xy_models, name, *stats)
        for name, stats in XY_MODEL_STATS.items()
    },
    "box-hull-dist-errors": check_box_hull_dist_errors,
    "cumulative-precision-recall": check_cumulative_precision_recall,
    "roc-models": check_roc_models,
    "rolling-mae-vs-hull-dist": check_rolling_mae,
    "hist-clf-pred-hull-dist": check_hist_clf,
    "element-prevalence-vs-error": check_element_prevalence,
    "scatter-largest-fp-diff-each-error": check_scatter_largest_fp_diff,
    "hist-largest-each-errors-fp-diff": check_hist_largest_each_errors,
    "hist-wbm-e-form-per-atom": check_hist_wbm_e_form_per_atom,
    "hist-wbm-hull-dist": check_hist_wbm_hull_dist,
    "spacegroup-sunbursts": check_spacegroup_sunbursts,
    "mp-vs-mp-trj-vs-wbm-arity-hist": check_arity_hist,
    "mp-trj-hists": check_mp_trj_hists,
    "mp-elemental-ref-energies": check_mp_elemental_ref_energies,
    "element-counts-mp-vs-wbm": check_element_counts,
    "spg-sankeys": check_spg_sankeys,
    "kappa-103-analysis": check_kappa_103_analysis,
}


@pytest.mark.parametrize("name", EXPECTED_PAYLOADS)
def test_payload_shape(name: str) -> None:
    """Each committed payload matches the shape its consumer page expects."""
    EXPECTED_PAYLOADS[name]()


# === staleness guards ===
# catch partial/stale payloads: a non-full regen silently drops models (as happened
# before the geo-opt CDF was refreshed from 7 to its full roster)

# all discovery figures share one roster: every active model with discovery metrics
DISCOVERY_PAYLOADS = (
    "box-hull-dist-errors",
    "cumulative-precision-recall",
    "roc-models",
    "rolling-mae-vs-hull-dist",
    "hist-clf-pred-hull-dist",
)


@pytest.mark.parametrize("name", DISCOVERY_PAYLOADS)
def test_discovery_payload_covers_active_models(name: str) -> None:
    """Each discovery figure covers every active model with discovery metrics."""
    expected = {model.key for model in Model.active() if model.metrics.get("discovery")}
    assert len(expected) > 30, f"sanity: too few discovery models ({len(expected)})"
    keys = payload_model_keys(name)
    assert keys == expected, (
        f"{name} roster drift: missing={expected - keys}, extra={keys - expected}. "
        "Run `uv run --with-editable . scripts/ingest_model.py <your-model> "
        "--payloads-only` to splice your model's entries into the committed "
        "payloads (needs only your own prediction files)."
    )


GEO_OPT_PAYLOADS = ("spg-sankeys", "struct-rmsd-cdf", "sym-ops-diff-bar")
TMI_PAYLOADS = (
    "hist-largest-each-errors-fp-diff",
    "scatter-largest-each-errors-fp-diff",
    "scatter-largest-fp-diff-each-error",
)
ELEMENT_PAYLOADS = {"element-prevalence-vs-error", "per-element-each-errors"}
PREDICTION_ROLES = {"generator", "payload_numerics"}
ERROR_ROLES = PREDICTION_ROLES | {"prediction_error_loader"}
RECIPE_ROLES = {
    "box-hull-dist-errors": PREDICTION_ROLES,
    "cumulative-precision-recall": PREDICTION_ROLES | {"stability_metrics"},
    "hist-clf-pred-hull-dist": PREDICTION_ROLES | {"stability_metrics"},
    "roc-models": PREDICTION_ROLES,
    "rolling-mae-vs-hull-dist": PREDICTION_ROLES,
    "spg-sankeys": {"generator", "payload_numerics"},
    "struct-rmsd-cdf": {"generator", "payload_numerics"},
    "sym-ops-diff-bar": {"generator"},
    **dict.fromkeys(TMI_PAYLOADS, ERROR_ROLES),
    "element-prevalence-vs-error": ERROR_ROLES | {"element_error_analysis"},
    "per-element-each-errors": {"generator", "prediction_error_loader"},
    "kappa-103-analysis": {
        "generator",
        "kappa_metrics",
        "payload_numerics",
        "phonon_reader",
        "phonon_schema",
        "thermal_conductivity",
    },
}


def test_multi_model_payload_provenance_matches_computation() -> None:
    """Every JSONL validates and hashes exactly its current computation sources."""
    paths = {
        name.removesuffix(".jsonl"): f"{SITE_FIG_DATA}/{name}"
        for name in os.listdir(SITE_FIG_DATA)
        if name.endswith(".jsonl")
    }
    paths["per-element-each-errors"] = (
        f"{SITE_DIR}/routes/models/per-element-each-errors.jsonl"
    )
    assert set(paths) == set(RECIPE_ROLES)
    for name, path in paths.items():
        payload = figs.read_jsonl_payload(path)
        identity = payload["identity"]
        audit = payload["audit"]
        sources = identity["recipe"]["sources"]
        for source in sources:
            with open(f"{ROOT}/{source['path']}", "rb") as file:
                source_bytes = file.read()
            assert source["size"] == len(source_bytes), source["path"]
            assert source["sha256"] == hashlib.sha256(source_bytes).hexdigest(), name

        expected_roles = RECIPE_ROLES[name]
        roles = {source["role"] for source in sources}
        assert roles == expected_roles, name
        if name in set(TMI_PAYLOADS) | ELEMENT_PAYLOADS:
            is_element_payload = name in ELEMENT_PAYLOADS
            benchmark_roles = {item["role"] for item in identity["benchmark_inputs"]}
            expected_benchmarks = {"wbm_summary"} | (
                {"mp_element_occurrences"} if is_element_payload else set()
            )
            assert benchmark_roles == expected_benchmarks, name
            assert ("pymatgen" in audit["runtime"]["packages"]) is is_element_payload


@pytest.mark.parametrize("name", GEO_OPT_PAYLOADS)
def test_geo_opt_payload_covers_active_models(name: str) -> None:
    """Each geo-opt payload includes every active model with an analysis file."""
    expected = {
        model.key
        for model in Model.active()
        if (geo_opt := model.metrics.get("geo_opt") or {}).get("pred_file")
        and (geo_opt.get("symprec=1e-5") or {}).get("analysis_file")
    }
    assert len(expected) > 30, f"sanity: too few geo-opt models ({len(expected)})"
    keys = payload_model_keys(name)
    assert keys == expected, (
        f"{name} roster drift: missing={expected - keys}, extra={keys - expected}. "
        "Run `python scripts/evals/geo_opt.py --auto-download`."
    )


def test_kappa_payload_covers_active_models() -> None:
    """kappa-103-analysis covers every active model with kappa_103 predictions."""
    expected = {
        model.key
        for model in Model.active()
        if ((model.metrics.get("phonons") or {}).get("kappa_103") or {}).get(
            "pred_file"
        )
    }
    assert len(expected) > 30, f"sanity: too few kappa models ({len(expected)})"
    keys = payload_model_keys("kappa-103-analysis")
    assert keys == expected, (
        f"kappa-103-analysis roster drift: missing={expected - keys}, "
        f"extra={keys - expected}. Run `uv run --with-editable . "
        "scripts/ingest_model.py <your-model> --payloads-only` to splice your "
        "model's entries into the committed payloads."
    )


# Sibling figures from a shared data source must agree on their stable-key roster.
PAYLOAD_FAMILIES = {
    "tmi-extras": (
        "element-prevalence-vs-error",
        "scatter-largest-fp-diff-each-error",
        "scatter-largest-each-errors-fp-diff",
        "hist-largest-each-errors-fp-diff",
    ),
}


@pytest.mark.parametrize("family", PAYLOAD_FAMILIES)
def test_payload_family_roster_consistent(family: str) -> None:
    """Payloads in a family share a model roster; a partial regen of one fails here."""
    names = PAYLOAD_FAMILIES[family]
    base = payload_model_keys(names[0])
    for name in names[1:]:
        roster = payload_model_keys(name)
        assert roster == base, f"{name} roster drifts from {names[0]}: {roster ^ base}"
