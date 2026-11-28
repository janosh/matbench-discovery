"""Tests for scripts/analyze_geo_opt.py's per-model, per-symprec analysis step."""

import gzip
import json
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Barrier
from types import ModuleType
from unittest.mock import PropertyMock, patch

import pandas as pd
import pytest
from pymatgen.core import Lattice, Structure

from matbench_discovery.enums import Key, Model
from tests.utils import import_repo_script

MOYO_VERSION = "0.4.2"
SYMPREC = 1e-2
ANALYSIS_CSV = f"2026-07-01-geo-opt-symprec=1e-2-moyo={MOYO_VERSION}.csv.gz"
# registered analysis from an older moyopy: must be recomputed, never overwritten
STALE_CSV = "2026-07-01-geo-opt-symprec=1e-2-moyo=0.4.1.csv.gz"


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
            record = {"material_id": mat_id, "structure": struct.as_dict()}
            file.write(f"{json.dumps(record)}\n")
    return geo_opt_path


def fake_analysis(material_ids: list[str]) -> pd.DataFrame:
    """Analysis frame with the columns calc_geo_opt_metrics needs."""
    return pd.DataFrame(
        {
            Key.structure_rmsd_vs_dft: [0.01] * len(material_ids),
            Key.spg_num_diff: [0] * len(material_ids),
            Key.n_sym_ops_diff: [0] * len(material_ids),
        },
        index=pd.Index(material_ids, name="material_id"),
    )


@pytest.mark.parametrize(
    ("debug_mode", "cached", "geo_opt_metrics"),
    [
        (0, True, {}),
        # a null symprec entry in the model YAML reads as no published analysis
        (0, True, {"symprec=1e-2": None}),
        (
            0,
            False,
            {"symprec=1e-2": {"analysis_file": {"name": f"model/{STALE_CSV}"}}},
        ),
        (3, False, {}),
        (3, True, {}),
    ],
    ids=[
        "cached-full",
        "cached-full-null-yaml",
        "stale-moyo-version",
        "debug-fresh",
        "debug-cached",
    ],
)
def test_analyze_model_symprec_metrics_and_debug_paths(
    debug_mode: int,
    cached: bool,
    geo_opt_metrics: dict[str, object],
    analyze_geo_opt: ModuleType,
    tmp_path: Path,
) -> None:
    """Cached CSVs skip the relaxed-structure file; debug runs spare real outputs."""
    material_ids = [f"wbm-1-{idx}" for idx in range(5)]
    geo_opt_path = write_geo_opt_file(tmp_path / "model", material_ids)
    old_path = tmp_path / "model" / STALE_CSV
    fake_analysis(material_ids[:1]).to_csv(old_path)
    historical_bytes = old_path.read_bytes()
    full_csv_path = tmp_path / "model" / ANALYSIS_CSV
    debug_csv_path = tmp_path / "model" / "debug" / str(debug_mode) / ANALYSIS_CSV
    expected_csv_path = debug_csv_path if debug_mode else full_csv_path
    n_analyzed = debug_mode or len(material_ids)
    if cached:
        expected_csv_path.parent.mkdir(parents=True, exist_ok=True)
        fake_analysis(material_ids[:n_analyzed]).to_csv(expected_csv_path)
    if debug_mode:
        # a cached CSV from a different subset size must not be picked up
        other_debug_csv = tmp_path / "model" / "debug" / "10" / ANALYSIS_CSV
        other_debug_csv.parent.mkdir(parents=True)
        fake_analysis(material_ids[:4]).to_csv(other_debug_csv)

    def fake_pred_vs_ref(
        df_model_analysis: pd.DataFrame, *_args: object, **_kwargs: object
    ) -> pd.DataFrame:
        return fake_analysis(list(df_model_analysis.index))

    def fake_sym_info(structs: dict[str, Structure], **_kwargs: object) -> pd.DataFrame:
        return pd.DataFrame(index=pd.Index(structs, name="material_id"))

    pred_file = {"name": os.path.relpath(geo_opt_path, tmp_path)}
    metrics = {"geo_opt": {"pred_file": pred_file, **geo_opt_metrics}}
    with (
        patch.object(analyze_geo_opt, "ROOT", str(tmp_path)),
        patch.object(
            Model,
            "geo_opt_path",
            new_callable=PropertyMock,
            return_value=str(geo_opt_path),
        ) as mock_geo_opt_path,
        patch.object(Model, "metrics", new_callable=PropertyMock, return_value=metrics),
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
            Model["alignn"],
            SYMPREC,
            MOYO_VERSION,
            df_dft_analysis=pd.DataFrame(index=pd.Index(material_ids)),
            dft_structs={},
            debug_mode=debug_mode,
        )

    assert df_metrics is not None
    assert old_path.read_bytes() == historical_bytes
    # metrics are derived from the (cached or fresh) analysis in every mode
    assert df_metrics["n_structures"].iloc[0] == n_analyzed
    # full runs score the whole test set, debug runs only their subset: 1 - 0.01/0.05
    assert df_metrics["rmsd_auc"].iloc[0] == pytest.approx(0.8, rel=1e-12)
    assert mock_sym_info.call_count == (0 if cached else 1)
    # a cached analysis must not resolve (download) or parse the multi-GB file
    assert mock_geo_opt_path.call_count == (0 if cached else 1)
    assert os.path.isfile(expected_csv_path)
    if debug_mode:
        # debug subsets land in debug/ and are never written to the model YAML
        assert not os.path.isfile(full_csv_path)
        mock_write.assert_not_called()
    else:
        mock_write.assert_called_once()
        assert mock_write.call_args.kwargs["regenerated"] is not cached
        # script joins with "/", so normalize before comparing (Windows uses "\")
        written_csv_path = os.path.normpath(mock_write.call_args.args[3])
        assert written_csv_path == os.path.normpath(full_csv_path)


def test_write_dft_reference_never_overwrites(
    analyze_geo_opt: ModuleType, tmp_path: Path
) -> None:
    """--write-dft-reference writes a new dated file, never clobbering an old one."""
    df_sym = pd.DataFrame(
        {"spg_num": [225]}, index=pd.Index(["wbm-1-1"], name="material_id")
    )
    barrier = Barrier(2)
    args = ({}, 1e-2, MOYO_VERSION, str(tmp_path))

    def simultaneous(*_args: object, **_kwargs: object) -> pd.DataFrame:
        """Let both workers pass the existing-file check before publication."""
        barrier.wait(timeout=10)
        return df_sym

    with (
        patch.object(
            analyze_geo_opt.symmetry, "get_sym_info_from_structs", simultaneous
        ),
        ThreadPoolExecutor(max_workers=2) as pool,
    ):
        jobs = [
            pool.submit(analyze_geo_opt.write_dft_reference, *args) for _ in range(2)
        ]
        errors = [job.exception() for job in jobs]
        assert sum(error is None for error in errors) == 1
        assert sum(isinstance(error, FileExistsError) for error in errors) == 1
        with pytest.raises(FileExistsError, match="Refusing to overwrite"):
            analyze_geo_opt.write_dft_reference(*args)
    (out_path,) = (tmp_path / "wbm").iterdir()
    assert out_path.name.endswith(
        f"-dft-geo-opt-symprec=1e-2-moyo={MOYO_VERSION}.csv.gz"
    )
    pd.testing.assert_frame_equal(
        pd.read_csv(out_path, index_col="material_id"), df_sym
    )


@pytest.mark.parametrize(
    ("moyo_version", "write_reference", "error"),
    [
        ("0.21.0", False, "WBM load reached"),
        ("0.22.0", False, "moyopy 0.22.0 does not match the registered"),
        # substring matches must not pass (0.2 is not 0.21.0)
        ("0.2", False, "moyopy 0.2 does not match the registered"),
        ("0.22.0", True, "WBM load reached"),
    ],
)
def test_main_checks_dft_reference_against_moyopy(
    moyo_version: str,
    write_reference: bool,
    error: str,
    analyze_geo_opt: ModuleType,
) -> None:
    """main() compares registered DFT references to moyopy before the WBM load."""
    raw_args = ["--models", "alignn", "--symprec", "1e-2"]
    raw_args += ["--write-dft-reference"] * write_reference
    with (
        patch.object(
            analyze_geo_opt.importlib.metadata, "version", return_value=moyo_version
        ),
        patch.object(
            analyze_geo_opt.DataFiles,
            "path",
            new_callable=PropertyMock,
            side_effect=RuntimeError("WBM load reached"),
        ),
        patch.object(
            analyze_geo_opt.DataFiles,
            "rel_path",
            new_callable=PropertyMock,
            return_value="wbm/dft-geo-opt-symprec=1e-2-moyo=0.21.0.csv.gz",
        ),
        pytest.raises((ValueError, RuntimeError), match=error),
    ):
        analyze_geo_opt.main(raw_args)
