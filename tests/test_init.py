from __future__ import annotations

import os

from matbench_discovery import ROOT, timestamp, today


def test_has_globals() -> None:
    assert os.path.isdir(ROOT)
    assert today == timestamp.split("@")[0]
    assert len(timestamp) == 19
