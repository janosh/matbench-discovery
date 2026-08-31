"""Tests for scripts/analyze_geo_opt.py's per-model, per-symprec analysis step."""

import gzip
import json
import os
from pathlib import Path
from types import ModuleType
from unittest.mock import PropertyMock, patch

import pandas as pd
import pytest
from pymatgen.core import Lattice, Structure

from matbench_discovery.enums import MbdKey, Model
from tests.utils import import_repo_script

MOYO_VERSION = "0.4.2"
SYMPREC = 1e-2


@pytest.fixture(scope="module", name="analyze_geo_opt")
def analyze_geo_opt_fixture() -> ModuleType:
    """scripts/analyze_geo_opt.py imported once for all tests in this module."""
    return import_repo_script("analyze_geo_opt", "scripts/analyze_geo_opt.py")


def write_geo_opt_file(model_dir: Path, material_ids: list[str]) -> Path:
    """Write a canonical geo-opt JSONL artifact with one cubic structure per ID."""
    struct = Structure(Lattice.cubic(3.0), ["Na", "Cl"], [[0] * 3, [0.5] * 3])
    geo_opt_path = model_dir / "2026-07-01-geo-opt.jsonl.gz"
    model_dir.mkdir(parents=True, exist_ok=True)
    with gzip.open(geo_opt_path, mode="wt") as file:
        for mat_id in material_ids:
            file.write(
                json.dumps({"material_id": mat_id, "structure": struct.as_dict()})
            )
            file.write("\n")
    return geo_opt_path


def fake_analysis(material_ids: list[str]) -> pd.DataFrame:
    """Analysis frame with the columns calc_geo_opt_metrics needs."""
    return pd.DataFrame(
        {
            MbdKey.structure_rmsd_vs_dft: [0.1] * len(material_ids),
            MbdKey.spg_num_diff: [0] * len(material_ids),
            MbdKey.n_sym_ops_diff: [0] * len(material_ids),
        },
        index=pd.Index(material_ids, name="material_id"),
    )


@pytest.mark.parametrize(
    ("debug_mode", "cached"),
    [(0, True), (3, False), (3, True)],
    ids=["cached-full", "debug-fresh", "debug-cached"],
)
def test_analyze_model_symprec_metrics_and_debug_paths(
    debug_mode: int,
    cached: bool,
    analyze_geo_opt: ModuleType,
    tmp_path: Path,
) -> None:
    """Cached CSVs still yield metrics; debug runs never touch the real outputs."""
    material_ids = [f"wbm-1-{idx}" for idx in range(5)]
    geo_opt_path = write_geo_opt_file(tmp_path / "model", material_ids)
    analysis_name = f"2026-07-01-geo-opt-symprec=1e-2-moyo={MOYO_VERSION}.csv.gz"
    full_csv_path = tmp_path / "model" / analysis_name
    debug_csv_path = tmp_path / "model" / "debug" / analysis_name
    expected_csv_path = debug_csv_path if debug_mode else full_csv_path
    n_analyzed = debug_mode or len(material_ids)
    if cached:
        expected_csv_path.parent.mkdir(parents=True, exist_ok=True)
        fake_analysis(material_ids[:n_analyzed]).to_csv(expected_csv_path)

    def fake_pred_vs_ref(
        df_model_analysis: pd.DataFrame, *_args: object, **_kwargs: object
    ) -> pd.DataFrame:
        return fake_analysis(list(df_model_analysis.index))

    def fake_sym_info(structs: dict[str, Structure], **_kwargs: object) -> pd.DataFrame:
        return pd.DataFrame(index=pd.Index(structs, name="material_id"))

    with (
        patch.object(
            Model,
            "geo_opt_path",
            new_callable=PropertyMock,
            return_value=str(geo_opt_path),
        ),
        patch.object(Model, "metrics", new_callable=PropertyMock, return_value={}),
        patch.object(
            analyze_geo_opt.symmetry,
            "get_sym_info_from_structs",
            side_effect=fake_sym_info,
        ) as mock_sym_info,
        patch.object(
            analyze_geo_opt.symmetry,
            "pred_vs_ref_struct_symmetry",
            side_effect=fake_pred_vs_ref,
        ),
        patch.object(analyze_geo_opt.geo_opt, "write_metrics_to_yaml") as mock_write,
    ):
        df_metrics = analyze_geo_opt.analyze_model_symprec(
            Model.alignn,
            SYMPREC,
            MOYO_VERSION,
            df_dft_analysis=pd.DataFrame(),
            dft_structs={},
            debug_mode=debug_mode,
        )

    assert df_metrics is not None
    # metrics are derived from the (cached or fresh) analysis in every mode
    assert df_metrics["n_structures"].iloc[0] == n_analyzed
    assert mock_sym_info.call_count == (0 if cached else 1)
    assert os.path.isfile(expected_csv_path)
    if debug_mode:
        # debug subsets land in debug/ and are never written to the model YAML
        assert not os.path.isfile(full_csv_path)
        mock_write.assert_not_called()
    else:
        mock_write.assert_called_once()
        assert mock_write.call_args.args[3] == str(full_csv_path)
