"""The shared calculator path must reach the specialist's precision guard."""

from pathlib import Path

import pytest

from matbench_discovery.calculators import CALCULATORS


def test_pet_pair_shared_registry_reaches_precision_guard(tmp_path: Path) -> None:
    """Missing registration or dropped dtype bypasses the frozen CPU/fp32 contract."""
    assert "pet_xs_pair_bo_v1" in CALCULATORS
    with pytest.raises(ValueError, match="float32"):
        CALCULATORS["pet_xs_pair_bo_v1"].make_calc(
            "cpu", checkpoint=str(tmp_path / "not-read.json"), dtype="float64"
        )
