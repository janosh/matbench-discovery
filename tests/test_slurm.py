import sys
from unittest.mock import patch

import pytest

from mb_discovery.slurm import _get_calling_file_path, slurm_submit_python


def test_slurm_submit() -> None:

    sys.argv += ["slurm-submit"]
    with pytest.raises(SystemExit) as exc_info, patch(
        "mb_discovery.slurm.subprocess.run"
    ) as mock_run:
        slurm_submit_python(
            job_name="test_job",
            log_dir="test_log_dir",
            time="0:0:1",
            slurm_flags=("--test_flag",),
        )
        assert exc_info.value.code == 0
        assert mock_run.call_count == 1


def test_get_calling_file_path() -> None:
    assert _get_calling_file_path(frame=1) == __file__

    def wrapper(frame: int) -> str:
        return _get_calling_file_path(frame)

    assert wrapper(frame=2) == __file__
