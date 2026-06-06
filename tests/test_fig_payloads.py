"""Shape guard for the committed data-only figure payloads in site/src/figs.

The Svelte pages cast these payloads to the interfaces in site/src/lib/fig-types.ts
at build time, so shape drift between the Python exporters and the site only
surfaces as broken figures. These tests mirror fig-types.ts and fail fast instead.
"""

from __future__ import annotations

import gzip
import json
import os
from functools import partial
from typing import Any

import pytest

from matbench_discovery import SITE_FIG_DATA


def load_payload(name: str) -> Any:
    """Load a committed figure payload by file stem."""
    with gzip.open(f"{SITE_FIG_DATA}/{name}.json.gz") as file:
        # our exporter writes allow_nan=False; reject NaN/Infinity literals that
        # python's json would otherwise happily parse back
        return json.load(
            file, parse_constant=lambda const: pytest.fail(f"non-finite {const=}")
        )


def assert_num_list(values: Any, *, length: int | None = None) -> None:
    """Assert a list of finite numbers (None allowed for gaps), optionally sized."""
    assert isinstance(values, list), f"expected list, got {type(values)=}"
    assert values, "expected non-empty list"
    assert all(val is None or isinstance(val, (int, float)) for val in values)
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
    for entry in models:
        assert isinstance(entry["label"], str)
        assert entry["label"]
        for key in keys:
            assert key in entry, f"model {entry['label']!r} missing {key!r}"
    return models


def test_no_orphan_payloads() -> None:
    """Every committed payload file is covered by a shape test below (and thus has a
    consumer page); orphans should be deleted, not committed.
    """
    on_disk = {
        file.removesuffix(".json.gz")
        for file in os.listdir(SITE_FIG_DATA)
        if file.endswith(".json.gz")
    }
    assert on_disk == set(EXPECTED_PAYLOADS), (
        f"unexpected={on_disk - set(EXPECTED_PAYLOADS)}, "
        f"missing={set(EXPECTED_PAYLOADS) - on_disk}"
    )


def check_box_hull_dist_errors() -> None:
    for model in assert_models(load_payload("box-hull-dist-errors"), "color"):
        assert_num_list(model["quantiles"], length=5)  # q05, q25, median, q75, q95


def check_cumulative_precision_recall() -> None:
    payload = load_payload("cumulative-precision-recall")
    assert payload["n_stable"] > 10_000
    for model in assert_models(payload, "color"):
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
    for model in assert_models(payload, "color"):
        assert_num_list(model["y"], length=len(payload["x"]))
    assert_xy(payload["density"])


def check_hist_clf() -> None:
    payload = load_payload("hist-clf-pred-hull-dist")
    assert_num_list(payload["bin_centers"])
    n_bins = len(payload["bin_centers"])
    for model in assert_models(payload, "f1"):
        for clf in ("tp", "fn", "fp", "tn"):
            assert_num_list(model[clf], length=n_bins)


def check_element_prevalence() -> None:
    payload = load_payload("element-prevalence-vs-error")
    elements = payload["elements"]
    assert all(isinstance(el, str) for el in elements)
    assert_num_list(payload["occurrences"], length=len(elements))
    for model in assert_models(payload, "color"):
        assert_num_list(model["y"], length=len(elements))


def check_scatter_largest_fp_diff() -> None:
    payload = load_payload("scatter-largest-fp-diff-each-error")
    assert_num_list(payload["fp_diff"])
    for model in assert_models(payload, "mae", "color"):
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
        roots = payload[key]
        assert len(roots) >= 6  # crystal systems
        for node in roots:
            assert node["label"]
            assert node["value"] > 0
            assert isinstance(node["children"], list)


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
        nodes, links = model["nodes"], model["links"]
        assert all(node["label"] for node in nodes)
        for link in links:
            assert 0 <= link["source"] < len(nodes)
            assert 0 <= link["target"] < len(nodes)
            assert link["value"] > 0


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
}


@pytest.mark.parametrize("name", EXPECTED_PAYLOADS)
def test_payload_shape(name: str) -> None:
    """Each committed payload matches the shape its consumer page expects."""
    EXPECTED_PAYLOADS[name]()
