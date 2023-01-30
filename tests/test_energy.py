from __future__ import annotations

import math
from typing import Any, Callable

import numpy as np
import pandas as pd
import pytest
from pymatgen.analysis.phase_diagram import PDEntry
from pymatgen.core import Lattice, Structure
from pymatgen.entries.computed_entries import (
    ComputedEntry,
    ComputedStructureEntry,
    Entry,
)
from pytest import approx
from sklearn.metrics import classification_report

from matbench_discovery.energy import (
    classify_stable,
    get_e_form_per_atom,
    get_elemental_ref_entries,
    stable_metrics,
)

dummy_struct = Structure(
    lattice=Lattice.cubic(5),
    species=("Fe", "O"),
    coords=((0, 0, 0), (0.5, 0.5, 0.5)),
)


@pytest.mark.parametrize(
    "constructor", [PDEntry, ComputedEntry, ComputedStructureEntry, lambda **x: x]
)
def test_get_e_form_per_atom(
    constructor: Callable[..., Entry | dict[str, Any]]
) -> None:
    """Test that the formation energy of a composition is computed correctly."""

    entry = {"composition": {"Fe": 1, "O": 1}, "energy": -2.5}
    elemental_ref_entries = {
        "Fe": {"composition": {"Fe": 1}, "energy": -1.0},
        "O": {"composition": {"O": 1}, "energy": -1.0},
    }
    if constructor == ComputedStructureEntry:
        entry["structure"] = dummy_struct
        entry.pop("composition")

    entry = constructor(**entry)

    # don't use ComputedStructureEntry for elemental ref entries, would need many
    # dummy structures
    if constructor == ComputedStructureEntry:
        constructor = ComputedEntry
    elemental_ref_entries = {
        k: constructor(**v) for k, v in elemental_ref_entries.items()
    }
    assert get_e_form_per_atom(entry, elemental_ref_entries) == -0.25


@pytest.mark.parametrize("constructor", [PDEntry, ComputedEntry, lambda **x: x])
@pytest.mark.parametrize("verbose", [True, False])
def test_get_elemental_ref_entries(
    constructor: Callable[..., Entry | dict[str, Any]], verbose: bool
) -> None:
    """Test that the elemental reference entries are correctly identified."""
    entries = [
        ("Fe1 O1", -2.5),
        ("Fe1", -1.0),
        ("Fe1", -2.0),
        ("O1", -1.0),
        ("O3", -2.0),
    ]
    elemental_ref_entries = get_elemental_ref_entries(
        [constructor(composition=comp, energy=energy) for comp, energy in entries],
        verbose=verbose,
    )
    if constructor.__name__ == "<lambda>":
        expected = {"Fe": PDEntry(*entries[2]), "O": PDEntry(*entries[3])}
    else:
        expected = {"Fe": constructor(*entries[2]), "O": constructor(*entries[3])}

    assert elemental_ref_entries == expected


@pytest.mark.parametrize(
    "stability_threshold, expected",
    [(-0.1, (6, 6, 11, 7)), (0, (8, 5, 10, 7)), (0.1, (10, 4, 10, 6))],
)
def test_classify_stable(
    stability_threshold: float, expected: tuple[int, int, int, int]
) -> None:
    np.random.seed(0)  # for reproducibility, makeDataFrame() uses np.random
    df = pd.util.testing.makeDataFrame()

    true_pos, false_neg, false_pos, true_neg = classify_stable(
        e_above_hull_true=df.A,
        e_above_hull_pred=df.B,
        stability_threshold=stability_threshold,
    )
    n_true_pos, n_false_neg, n_false_pos, n_true_neg = map(
        sum, (true_pos, false_neg, false_pos, true_neg)
    )

    assert (n_true_pos, n_false_neg, n_false_pos, n_true_neg) == expected
    assert n_true_pos + n_false_neg + n_false_pos + n_true_neg == len(df)
    assert n_true_neg + n_false_pos == sum(df.A > stability_threshold)
    assert n_true_pos + n_false_neg == sum(df.A <= stability_threshold)


def test_stable_metrics() -> None:
    metrics = stable_metrics(np.arange(-1, 1, 0.1), np.arange(1, -1, -0.1))
    for key, val in dict(
        DAF=0,
        Precision=0,
        Recall=0,
        Prevalence=0.55,
        Accuracy=0,
        TPR=0,
        FPR=1,
        TNR=0,
        FNR=1,
        MAE=0.999,
        RMSE=1.157,
        R2=-3.030,
    ).items():
        assert metrics[key] == approx(val, abs=1e-3), f"{key=}"

    assert math.isnan(metrics["F1"])

    # test that feeding in random numpy data gives the same result as feeding it into
    # sklearn.metrics.classification_report
    np.random.seed(0)
    y_true, y_pred = np.random.randn(100, 2).T
    metrics = stable_metrics(y_true, y_pred)

    report = classification_report(y_true > 0, y_pred > 0, output_dict=True)
    assert metrics["Precision"] == report["False"]["precision"]
    assert metrics["Recall"] == report["False"]["recall"]
    assert metrics["F1"] == report["False"]["f1-score"]

    assert stable_metrics.__doc__  # for mypy
    assert all(key in stable_metrics.__doc__ for key in metrics)
