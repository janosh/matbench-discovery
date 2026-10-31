"""Tests for rerun scoring and verdict blocks in metrics.reproducibility."""

import gzip
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
import yaml
from jsonschema import Draft7Validator

from matbench_discovery import ROOT
from matbench_discovery.enums import Model
from matbench_discovery.metrics import reproducibility as repro
from tests.utils import patch_model_metadata

DIAT = {"a": 1.0, "b": 2.0}  # every diatomic metric decides the verdict
MD = {"rdf": 10.0, "energy_rmse": 3.0, "force_rmse": 40.0}
RUNS = [{"label": "rerun", "file": {"name": "a"}}]
with open(f"{ROOT}/tests/model-schema.yml", encoding="utf-8") as schema_file:
    DEFS: dict[str, Any] = yaml.safe_load(schema_file)["definitions"]
VALIDATOR = Draft7Validator(DEFS["Reproducibility"] | {"definitions": DEFS})


@pytest.mark.parametrize(
    ("task", "registered", "runs", "status", "rel_changes"),
    [
        ("diatomics", DIAT, [DIAT], "reproduced", {}),
        # max over runs at 3 sig figs, within tolerance; zero changes omitted
        ("diatomics", DIAT, [DIAT | {"a": 1.0005}, DIAT | {"a": 0.99921}],
         "reproduced", {"a": 7.9e-4}),
        ("diatomics", DIAT, [DIAT | {"b": 2.01}], "drifted", {"b": 5e-3}),
        # a nonzero rerun of a zero registered score is a 100% change
        ("diatomics", {"a": 0.0}, [{"a": 1.0}], "drifted", {"a": 1.0}),
        ("diatomics", DIAT, [{"a": 1.0}], "inconclusive", {}),
        # chaotic trajectory metrics are listed but never decide the verdict
        ("md", MD, [MD | {"rdf": 10.5, "energy_rmse": 3.00003}], "reproduced",
         {"rdf": 0.05, "energy_rmse": 1e-5}),
        ("md", MD, [MD | {"force_rmse": 40.2}], "drifted", {"force_rmse": 5e-3}),
        # frame-wise diagnostics missing from the registered file or a run
        ("md", {"rdf": 10.0}, [MD | {"rdf": 10.2}], "inconclusive", {"rdf": 0.02}),
        ("md", MD, [MD, {"rdf": 11.0}], "inconclusive", {"rdf": 0.1}),
    ],
)  # fmt: skip
def test_make_block_verdicts(
    task: repro.ReproTask,
    registered: dict[str, float],
    runs: list[dict[str, float]],
    status: str,
    rel_changes: dict[str, float],
) -> None:
    """Verdicts follow the per-task rules and blocks validate against the schema."""
    cause = "unknown" if status == "drifted" else None
    note = None if status == "reproduced" else "why"
    block = repro.make_block(task, registered, runs, runs=RUNS, cause=cause, note=note)
    assert block["status"] == status
    assert block["rel_changes"] == pytest.approx(rel_changes, rel=1e-12)
    assert list(block)[:2] == ["checked", "status"]
    assert not list(VALIDATOR.iter_errors(block))


@pytest.mark.parametrize(
    ("run", "status", "cause", "note"),
    [
        (DIAT, "reproduced", "unknown", None),
        (DIAT | {"a": 2.0}, "drifted", None, "why"),
        (DIAT | {"a": 2.0}, "drifted", "unknown", None),
        ({}, "inconclusive", None, None),
    ],
)
def test_make_block_cause_note_rules(
    run: dict[str, float], status: str, cause: str | None, note: str | None
) -> None:
    """Python and schema agree: a cause iff drifted, a note unless reproduced."""
    with pytest.raises(ValueError, match=f"{status} verdict needs"):
        repro.make_block("diatomics", DIAT, [run], runs=RUNS, cause=cause, note=note)
    block = repro.make_block("diatomics", DIAT, [DIAT], runs=RUNS)
    overrides = {"status": status, "cause": cause, "note": note}
    block |= {key: val for key, val in overrides.items() if val}
    assert list(VALIDATOR.iter_errors(block))


def test_make_block_repeat_rel_change() -> None:
    """Repeats measure the max change over ordered pairs of verdict metrics only."""
    repeats = [MD | {"rdf": 20.0, "force_rmse": 40.004}, MD]
    block = repro.make_block("md", MD, [MD], runs=RUNS, repeat_scores=repeats)
    assert block["repeat_rel_change"] == pytest.approx(1e-4, rel=1e-6)
    # 1.5 vs 1.0 -> 0.5 beats 1.0 vs 1.5 -> 1/3
    repeats = [DIAT | {"a": 1.2}, DIAT | {"a": 1.5}, DIAT]
    block = repro.make_block(
        "diatomics", DIAT, [DIAT], runs=RUNS, repeat_scores=repeats
    )
    assert block["repeat_rel_change"] == 0.5


def test_score_files(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Diatomics drop exclusions; MD drops cost, scales RMSEs and checks coverage."""
    model = Model.mace_mp_0
    dists = np.linspace(0.3, 6.0, 120)
    curves = {}
    for elem, depth in (("H", 4.5), ("Li", 1.0)):
        energies = depth * ((1 - np.exp(-1.5 * (dists - 1.0))) ** 2 - 1)
        forces = np.zeros((len(dists), 2, 3))
        forces[:, 0, 0] = -np.gradient(energies, dists)
        forces[:, 1, 0] = -forces[:, 0, 0]
        curves[elem] = {"energies": energies.tolist(), "forces": forces.tolist()}
    paths = [f"{tmp_path}/both.json.gz", f"{tmp_path}/h_only.json.gz"]
    for path, elems in zip(paths, (("H", "Li"), ("H",)), strict=True):
        with gzip.open(path, mode="wt", encoding="utf-8") as file:
            homo_nuclear = {elem: curves[elem] for elem in elems}
            json.dump({"distances": dists.tolist(), "homo-nuclear": homo_nuclear}, file)
    both, h_only = repro.score_files("diatomics", model, paths)
    assert both != h_only
    excluded = {"metrics": {"diatomics": {"excluded_formula_reasons": {"Li-Li": "x"}}}}
    patch_model_metadata(monkeypatch, model, excluded)
    assert repro.score_files("diatomics", model, paths[:1]) == [h_only]

    rows = {"rdf_error": [1.0, 3.0], "energy_rmse": [1e-3, 3e-3], "hardware": "H200"}
    df_md = pd.DataFrame(rows, index=pd.Index(["sys_a", "sys_b"], name="system"))
    df_md.to_csv(md_path := f"{tmp_path}/md.csv.gz")
    df_md.iloc[:1].to_csv(partial_path := f"{tmp_path}/partial.csv.gz")
    expected = {"rdf_error": 2.0, "energy_rmse": 2.0}
    assert repro.score_files("md", model, [md_path]) == [pytest.approx(expected)]
    with pytest.raises(ValueError, match=r"covers \['sys_a'\], expected"):
        repro.score_files("md", model, [md_path, partial_path])
