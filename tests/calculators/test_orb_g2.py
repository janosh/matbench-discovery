"""The standard registry stays dependency-free; optional runtime checks real guard."""

from pathlib import Path

import pytest

from matbench_discovery.calculators import CALCULATORS


def test_orb_g2_registration_without_runtime() -> None:
    assert "orb_g2_groupcv_v1" in CALCULATORS


def test_orb_g2_registry_reaches_precision_guard(tmp_path: Path) -> None:
    pytest.importorskip("orb_g2_cv", reason="optional model-specific runtime")
    with pytest.raises(ValueError, match="float32"):
        CALCULATORS["orb_g2_groupcv_v1"].make_calc(
            "cpu", checkpoint=str(tmp_path / "absent.json"), dtype="float64"
        )
