"""Shape guard for the committed data-only figure payloads in site/src/figs.

The Svelte pages type these payloads via the ambient module declarations in
site/src/figs/payloads.d.ts, so shape drift between the Python exporters and the site
only surfaces as broken figures. These tests mirror payloads.d.ts and fail fast.
"""

from __future__ import annotations

import gzip
import json
import math
import os
from functools import partial
from typing import Any

import pandas as pd
import pytest

from matbench_discovery import SITE_DIR, SITE_FIG_DATA, figs
from matbench_discovery.enums import Model
from scripts.model_figs.kappa_103_analysis import row_flag


def reject_json_constant(const: str) -> None:
    """Reject NaN/Infinity literals in JSON figure payloads."""
    raise ValueError(f"non-finite JSON constant {const!r}")


def load_payload(name: str) -> dict[str, Any]:
    """Load a committed figure payload by stem - either a gzipped aggregate
    ``<name>.json.gz`` or a line-delimited ``<name>.jsonl`` (reassembled like the site's
    jsonl Vite plugin). The exporter writes ``allow_nan=False`` so NaN can't reach disk;
    the aggregate path also rejects NaN/Infinity literals to catch drift.
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


def assert_models(payload: dict[str, Any], *keys: str, n_min: int = 2) -> list:
    """Assert a {models: [...]} payload where each entry has (at least) ``keys``."""
    models = payload["models"]
    assert len(models) >= n_min, f"expected >= {n_min} models, got {len(models)}"
    # ids must be unique - a git merge that duplicated a model line (or a generator bug)
    # would otherwise render two entries for one model
    ids = [entry.get("key") or entry["label"] for entry in models]
    assert len(ids) == len(set(ids)), f"duplicate model ids: {sorted(ids)}"
    for entry in models:
        assert isinstance(entry["label"], str)
        assert entry["label"]
        for key in keys:
            assert key in entry, f"model {entry['label']!r} missing {key!r}"
    return models


def payload_model_ids(name: str, field: str) -> set[str]:
    """Set of per-model identifiers (``field`` = 'key' or 'label') in a payload."""
    return {model[field] for model in load_payload(name)["models"]}


def test_no_orphan_payloads() -> None:
    """Every committed payload file is covered by a shape test below (and thus has a
    consumer page); orphans should be deleted, not committed.
    """
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
    """The route-local per-element-errors .jsonl (imported by per-element-errors.ts) is
    readable and well-shaped: each column maps a model_key/metadata label to a dict of
    finite per-element values (None allowed for gaps).
    """
    path = f"{SITE_DIR}/routes/models/per-element-each-errors.jsonl"
    columns = figs.read_jsonl_payload(path)["models"]
    assert len(columns) > 10, f"expected >10 columns, got {len(columns)}"
    for column in columns:
        assert isinstance(column["key"], str)
        assert_num_list(list(column["values"].values()))


def check_box_hull_dist_errors() -> None:
    for model in assert_models(load_payload("box-hull-dist-errors"), "key"):
        assert_num_list(model["quantiles"], length=5)  # q05, q25, median, q75, q95


def check_cumulative_precision_recall() -> None:
    payload = load_payload("cumulative-precision-recall")
    assert payload["n_stable"] > 10_000
    for model in assert_models(payload, "key"):
        assert_num_list(model["x"])
        assert_num_list(model["precision"], length=len(model["x"]))
        assert_num_list(model["recall"], length=len(model["x"]))
        assert_num_list(model["end"], length=3)


def check_roc_models() -> None:
    for model in assert_models(load_payload("roc-models"), "key", "auc"):
        assert 0.5 < model["auc"] <= 1
        assert_num_list(model["fpr"])
        assert_num_list(model["tpr"], length=len(model["fpr"]))


def check_rolling_mae() -> None:
    payload = load_payload("rolling-mae-vs-hull-dist")
    assert_num_list(payload["x"])
    for model in assert_models(payload, "key"):
        assert_num_list(model["y"], length=len(payload["x"]))
    assert_xy(payload["density"])


def check_hist_clf() -> None:
    payload = load_payload("hist-clf-pred-hull-dist")
    assert_num_list(payload["bin_centers"])
    n_bins = len(payload["bin_centers"])
    for model in assert_models(payload, "key", "f1"):
        for clf in ("tp", "fn", "fp", "tn"):
            assert_num_list(model[clf], length=n_bins)


def check_element_prevalence() -> None:
    payload = load_payload("element-prevalence-vs-error")
    elements = payload["elements"]
    assert all(isinstance(el, str) for el in elements)
    assert_num_list(payload["occurrences"], length=len(elements))
    for model in assert_models(payload):
        assert_num_list(model["y"], length=len(elements))


def check_scatter_largest_fp_diff() -> None:
    payload = load_payload("scatter-largest-fp-diff-each-error")
    assert_num_list(payload["fp_diff"])
    for model in assert_models(payload, "mae"):
        assert_num_list(model["y"], length=len(payload["fp_diff"]))


def check_hist_largest_each_errors() -> None:
    for model in assert_models(load_payload("hist-largest-each-errors-fp-diff")):
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
    for model in assert_models(load_payload("spg-sankeys"), "key"):
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
    for model in assert_models(payload, "key", "freq_w1_mean", "freq_pairs"):
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
    """Each discovery figure must include every active model with discovery metrics
    (by model_key), so a partial regen that drops models fails fast instead of silently
    shipping an incomplete leaderboard figure.
    """
    expected = {
        model.key
        for model in Model.active()
        if (disc := (model.metrics or {}).get("discovery")) and disc != "not applicable"
    }
    assert len(expected) > 30, f"sanity: too few discovery models ({len(expected)})"
    keys = payload_model_ids(name, "key")
    assert keys == expected, (
        f"{name} roster drift: missing={expected - keys}, extra={keys - expected}. "
        "Run `uv run --with-editable . scripts/ingest_model.py <your-model> "
        "--payloads-only` to splice your model's entries into the committed "
        "payloads (needs only your own prediction files)."
    )


def test_kappa_payload_covers_active_models() -> None:
    """The kappa-103-analysis payload must include every active model with kappa_103
    predictions, so a partial regen fails fast instead of silently dropping models.
    """
    expected = {
        model.key
        for model in Model.active()
        if isinstance(phonons := (model.metrics or {}).get("phonons"), dict)
        and (phonons.get("kappa_103") or {}).get("pred_file")
    }
    assert len(expected) > 30, f"sanity: too few kappa models ({len(expected)})"
    keys = payload_model_ids("kappa-103-analysis", "key")
    assert keys == expected, (
        f"kappa-103-analysis roster drift: missing={expected - keys}, "
        f"extra={keys - expected}. Run `uv run --with-editable . "
        "scripts/ingest_model.py <your-model> --payloads-only` to splice your "
        "model's entries into the committed payloads."
    )


# sibling figures from a shared data source must agree on their model roster (these
# carry only `label`, not `key`), so a partial regen of one fails against the rest
LABEL_PAYLOAD_FAMILIES = {
    "geo-opt": ("struct-rmsd-cdf", "sym-ops-diff-bar"),
    "tmi-extras": (
        "element-prevalence-vs-error",
        "scatter-largest-fp-diff-each-error",
        "scatter-largest-each-errors-fp-diff",
        "hist-largest-each-errors-fp-diff",
    ),
}


@pytest.mark.parametrize("family", LABEL_PAYLOAD_FAMILIES)
def test_label_payload_family_roster_consistent(family: str) -> None:
    """Payloads in a family share a model roster; a partial regen of one fails here."""
    names = LABEL_PAYLOAD_FAMILIES[family]
    base = payload_model_ids(names[0], "label")
    for name in names[1:]:
        roster = payload_model_ids(name, "label")
        assert roster == base, f"{name} roster drifts from {names[0]}: {roster ^ base}"
