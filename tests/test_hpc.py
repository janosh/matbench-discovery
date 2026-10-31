"""Tests for high-performance computing utilities."""

import json
import os
import sys
from functools import partial
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import pytest

from matbench_discovery import hpc


@pytest.mark.parametrize(
    (
        "py_file_path",
        "time",
        "account",
        "partition",
        "array",
        "pre_cmd",
        "submit_as_temp_file",
    ),
    [
        (None, None, None, None, None, "", True),
        (
            "path/to/file.py",
            "0:0:1",
            "fake-account",
            "fake-partition",
            None,
            "module load pytorch;",
            False,
        ),
        (None, "0:0:1", None, None, "0-9%2", "ENV_VAR=42", True),
        ("path/to/file.py", None, "fake-account", None, "1-10", "", False),
        (None, None, None, "fake-partition", None, "ENV_VAR=42", False),
    ],
)
def test_slurm_submit(
    capsys: pytest.CaptureFixture[str],
    tmp_path: Path,
    py_file_path: str | None,
    time: str | None,
    account: str | None,
    partition: str | None,
    array: str | None,
    pre_cmd: str,
    submit_as_temp_file: bool,
) -> None:
    """slurm_submit is a no-op outside submit mode and builds the sbatch command."""
    job_name, out_dir = "test_job", f"{tmp_path}/slurm-logs"
    kwargs = dict(
        job_name=job_name,
        out_dir=out_dir,
        time=time,
        partition=partition,
        account=account,
        py_file_path=py_file_path,
        slurm_flags="--foo",
        array=array,
        pre_cmd=pre_cmd,
        submit_as_temp_file=submit_as_temp_file,
    )

    hpc.slurm_submit(**kwargs)  # ty: ignore[invalid-argument-type]
    assert capsys.readouterr() == ("", "")  # does nothing in normal mode
    assert os.path.isdir(out_dir)  # sbatch fails if the log dir is missing

    with patch.dict(os.environ, {"SLURM_JOB_ID": "1234"}, clear=True):
        slurm_vars = hpc.slurm_submit(**kwargs)  # ty: ignore[invalid-argument-type]
    expected_pre_cmd = (
        pre_cmd if pre_cmd.endswith(";") or not pre_cmd else f"{pre_cmd};"
    )
    assert slurm_vars == {
        "slurm_job_id": "1234",
        "slurm_flags": "--foo",
        **({"slurm_timelimit": time} if time else {}),
        **({"pre_cmd": expected_pre_cmd} if pre_cmd else {}),
    }

    # submit mode prints the cmd and calls subprocess.run()
    with (
        pytest.raises(SystemExit),
        patch("sys.argv", ["slurm-submit"]),
        patch("matbench_discovery.hpc.subprocess.run") as mock_subprocess_run,
        patch(
            "matbench_discovery.hpc.tempfile.mkdtemp",
            return_value="/tmp/slurm_job_123",
        ),
        patch("matbench_discovery.hpc.shutil.copy2") as mock_copy2,
        patch.dict(os.environ, {"SLURM_JOB_ID": "1234"}, clear=True),
    ):
        hpc.slurm_submit(**kwargs)  # ty: ignore[invalid-argument-type]
    assert mock_subprocess_run.call_count == 1

    expected_py_file_path = py_file_path or __file__
    if submit_as_temp_file:
        temp_path = f"/tmp/slurm_job_123/{os.path.basename(expected_py_file_path)}"
        mock_copy2.assert_called_once_with(expected_py_file_path, temp_path)
        expected_py_file_path = temp_path
    else:
        mock_copy2.assert_not_called()

    log_name = "slurm-%A-%a.log" if array else "slurm-%A.log"
    wrap_cmd = f"{expected_pre_cmd} python {expected_py_file_path}".strip()
    sbatch_cmd = (
        f"sbatch --job-name {job_name} --output {out_dir}/{log_name} --foo "
        f"--wrap {wrap_cmd}"
    ).replace(" --", "\n  --")
    for key, val in (
        ("time", time),
        ("account", account),
        ("partition", partition),
        ("array", array),
    ):
        if val is not None:
            sbatch_cmd += f"\n  --{key} {val}"
    stdout, stderr = capsys.readouterr()
    assert sbatch_cmd in stdout
    assert stderr == ""


def test_get_calling_file_path() -> None:
    assert hpc._get_calling_file_path(frame=1) == __file__  # noqa: SLF001

    def wrapper(frame: int) -> str:
        return hpc._get_calling_file_path(frame)  # noqa: SLF001

    assert wrapper(frame=2) == __file__


@pytest.mark.parametrize(
    ("units", "expected"),
    [
        (
            [
                # unrelated keys like excluded_formula_reasons are ignored by merge
                {
                    "hardware": "H200",
                    "run_time_sec": 100,
                    "excluded_formula_reasons": {},
                },
                {"hardware": "H200", "run_time_sec": 50.5},
            ],
            {"hardware": "H200", "run_time_sec": 150.5},
        ),
        ([], {}),
        ([{}, {}], {}),  # no provenance -> empty (leaves existing YAML untouched)
        # all-or-nothing: a field missing from any unit is dropped, not summed
        # over the units that have it
        (
            [
                {"hardware": "H200", "run_time_sec": 10.0, "max_rss_gb": 4.0},
                {"hardware": "H200", "run_time_sec": 20.0},
            ],
            {"hardware": "H200", "run_time_sec": 30.0},
        ),
        (
            # times, loads and process counts sum (a process charges its load to
            # one unit), memory peaks max
            [
                {"hardware": "H200", "run_time_sec": 10.0, "load_time_sec": 3.0}
                | {"n_processes": 1, "max_rss_gb": 4.2, "max_gpu_mem_gb": 9.0},
                {"hardware": "H200", "run_time_sec": 20.0, "load_time_sec": 0.0}
                | {"n_processes": 0, "max_rss_gb": 6.1, "max_gpu_mem_gb": 11.5},
            ],
            {"hardware": "H200", "run_time_sec": 30.0, "load_time_sec": 3.0}
            | {"n_processes": 1, "max_rss_gb": 6.1, "max_gpu_mem_gb": 11.5},
        ),
        (
            # GPU-dependent costs need a hardware label: unlabeled units keep only
            # host RSS and process counts
            [
                {"run_time_sec": 10.0, "max_rss_gb": 4.2, "max_gpu_mem_gb": 11.5},
                {"run_time_sec": 20.0, "max_rss_gb": 6.1, "max_gpu_mem_gb": 9.0},
            ],
            {"max_rss_gb": 6.1},
        ),
        # NaN counts as missing, so the result can't depend on unit order
        ([{"max_rss_gb": float("nan")}, {"max_rss_gb": 2.0}], {}),
        ([{"max_rss_gb": 2.0}, {"max_rss_gb": float("nan")}], {}),
        # versions and labels only when every unit agrees (lists by value)
        ([{"versions": {"numpy": "1"}}, {}], {}),
        (
            [{"versions": {"numpy": "1"}, "acceleration": ["cuequivariance"]}] * 2,
            {"versions": {"numpy": "1"}, "acceleration": ["cuequivariance"]},
        ),
        (
            # latest completion time and sorted, de-duplicated host/job unions
            [
                {
                    "completed_at": "2026-07-02T00:00:00",
                    "hostname": "node-b",
                    "slurm_job_id": "2",
                    "slurm_array_task_id": "1",
                },
                {
                    "completed_at": "2026-07-03T00:00:00",
                    "hostname": "node-a",
                    "slurm_job_id": "2",
                    "slurm_array_task_id": "0",
                },
                {"completed_at": "2026-07-01T00:00:00", "hostname": "node-b"},
            ],
            {
                "completed_at": "2026-07-03T00:00:00",
                "hostnames": ["node-a", "node-b"],
                "slurm_job_ids": ["2"],
                "slurm_array_task_ids": ["0", "1"],
            },
        ),
    ],
)
def test_merge_audit_metadata(
    units: list[dict[str, object]], expected: dict[str, object]
) -> None:
    """Cost fields merge all-or-nothing; audit fields take latest time + unions."""
    assert hpc.merge_audit_metadata(units) == expected


@pytest.mark.parametrize(
    "changed",
    [
        {"hardware": ""},
        {"hardware": "A100"},
        {"acceleration": []},
        {"cost_source": None},
        {"n_shards": 2},
    ],
)
def test_merge_audit_metadata_mixed_labels_drop_cost(
    changed: dict[str, object],
) -> None:
    """Mixed run labels drop the cost (warning); mixed versions raise if strict."""
    labels = {"cost_source": "runner", "hardware": "H200", "acceleration": ["cueq"]}
    unit = labels | {"n_shards": 1, "run_time_sec": 1.0, "versions": {"numpy": "1"}}
    assert hpc.merge_audit_metadata([unit] * 2, strict=True) == unit | {
        "run_time_sec": 2.0
    }
    with pytest.warns(UserWarning, match="Dropping the whole cost block"):
        merged = hpc.merge_audit_metadata([unit, unit | changed], strict=True)
    assert merged == {"versions": {"numpy": "1"}}
    with pytest.raises(ValueError, match="mixed or missing package versions"):
        hpc.merge_audit_metadata(
            [unit, unit | {"versions": {"numpy": "2"}}], strict=True
        )


# flat artifact run metadata and the YAML cost mapping it yields
RUNNER_RUN = {"cost_source": "runner", "hardware": "NVIDIA H200", "run_time_sec": 5.0}
RUNNER_COST = {"source": "runner", "hardware": "NVIDIA H200", "run_time_sec": 5.0}
AUTHOR_COST = {"source": "author_reported", "hardware": "NVIDIA A100"}


@pytest.mark.parametrize(
    ("existing", "artifact", "same_artifact", "replace", "expected"),
    [
        # artifact cost wins whole: no stale key (max_rss_gb) survives next to it
        ({**RUNNER_COST, "max_rss_gb": 9.0}, RUNNER_RUN, True, False, RUNNER_COST),
        # legacy artifact (no cost_source) inherits an artifact-backed YAML tag and
        # contributes only its trustworthy keys (older runners recorded the
        # requested, not the actual, dtype)
        (
            {"source": "salvaged", "hardware": "NVIDIA H200"},
            {"hardware": "NVIDIA H200", "run_time_sec": 7.0, "dtype": "float64"},
            True,
            False,
            {"source": "salvaged", "hardware": "NVIDIA H200", "run_time_sec": 7.0},
        ),
        # author_reported blocks are kept untouched when recomputing an artifact
        # without cost (YAML-only author report) ...
        (AUTHOR_COST, {}, True, False, AUTHOR_COST),
        # ... re-derived but kept author_reported for an author-run artifact ...
        (
            AUTHOR_COST,
            {"hardware": "NVIDIA A100 80GB PCIe", "run_time_sec": 7.0},
            True,
            False,
            {**AUTHOR_COST, "hardware": "NVIDIA A100 80GB PCIe", "run_time_sec": 7.0},
        ),
        # ... replaced by a new run's cost only on explicit request ...
        (AUTHOR_COST, RUNNER_RUN, False, True, RUNNER_COST),
        (AUTHOR_COST, RUNNER_RUN, False, False, "refusing to replace author_reported"),
        # ... and dropped (warning) for a new artifact without cost
        (AUTHOR_COST, {}, False, False, "refusing to replace author_reported"),
        (AUTHOR_COST, {}, False, True, "dropping author_reported cost"),
        # a runner-tagged artifact overrides a salvaged tag on recompute
        ({"source": "salvaged", "hardware": "x"}, RUNNER_RUN, True, False, RUNNER_COST),
        # merges that dropped hardware/run time leave nothing publishable
        (RUNNER_COST, {"cost_source": "runner", "max_rss_gb": 1.0}, False, False, None),
        # untagged, contradicted or invalid cost raises
        (RUNNER_COST, {}, True, False, "records no cost"),
        ({"hardware": "NVIDIA H200"}, RUNNER_RUN, True, False, "no valid source"),
        ({}, {"hardware": "x", "run_time_sec": 1.0}, False, False, "no cost_source"),
        ({}, RUNNER_RUN | {"cost_source": "estimated"}, False, False, "must be one"),
    ],
)
def test_with_cost_block(
    existing: dict[str, object],
    artifact: dict[str, object],
    same_artifact: bool,
    replace: bool,
    expected: dict[str, object] | str | None,
) -> None:
    """Cost comes whole from artifact or YAML, never mixed or silently replaced."""
    section = {"F1": 0.5, "pred_file": "a"} | ({"cost": existing} if existing else {})
    write = partial(
        hpc.with_cost_block,
        section,
        artifact,
        after="pred_file",
        same_artifact=same_artifact,
        context="mace metrics.md",
        replace_author_reported=replace,
    )
    if isinstance(expected, str) and expected.startswith("dropping"):
        with pytest.warns(UserWarning, match=expected):
            assert "cost" not in write()
        return
    if isinstance(expected, str):
        with pytest.raises(ValueError, match=expected):
            write()
        return
    written = write()
    assert list(written) == ["F1", "pred_file", *(["cost"] if expected else [])]
    assert written.get("cost") == expected
    assert list(written.get("cost") or ()) == list(expected or ())  # canonical order


def test_detect_runtime_probe_fallback(monkeypatch: pytest.MonkeyPatch) -> None:
    """Calculators without reachable torch weights get the requested device."""
    for framework in ("jax", "tensorflow"):
        monkeypatch.delitem(sys.modules, framework, raising=False)
    monkeypatch.setattr(hpc, "detect_hardware", lambda _device: "NVIDIA H200")
    assert hpc.detect_runtime(object(), requested_device="cuda") == {
        "backend": "torch_probe",
        "device": "cuda",
        "torch_device": "cuda",
        "hardware": "NVIDIA H200",
    }


@pytest.mark.parametrize(
    ("runtime", "expected_gpu_gb", "succeeded_raises"),
    [
        # non-torch backends never report torch's (meaningless) CUDA counter
        ({"backend": "jax", "device": "gpu"}, None, False),
        ({"backend": "torch", "device": "cpu"}, None, False),
        ({"backend": "torch", "device": "cuda", "torch_device": "cuda:1"}, 2.5, False),
        # a failed unit may not have touched the GPU, but a torch GPU label
        # without VRAM after successful work raises
        ({"backend": "torch_probe", "device": "cuda"}, 0.0, True),
        ({"backend": "torch", "device": "cuda"}, 0.0, True),
    ],
)
def test_memory_peaks(
    runtime: dict[str, str], expected_gpu_gb: float | None, succeeded_raises: bool
) -> None:
    """memory_peaks reports RSS + weight-device VRAM, rejecting GPU labels w/o VRAM."""
    allocated = Mock(return_value=(expected_gpu_gb or 0) * 1e9)
    cuda = SimpleNamespace(is_available=lambda: True, max_memory_allocated=allocated)
    with patch.dict(sys.modules, torch=SimpleNamespace(cuda=cuda)):
        mem = hpc.memory_peaks(runtime, succeeded=False)
        if succeeded_raises:
            with pytest.raises(ValueError, match="allocated no CUDA memory"):
                hpc.memory_peaks(runtime, succeeded=True)
        else:
            assert hpc.memory_peaks(runtime, succeeded=True) == mem
    assert mem.get("max_gpu_mem_gb") == expected_gpu_gb
    if expected_gpu_gb is not None:
        allocated.assert_called_with(runtime.get("torch_device"))
    if sys.platform != "win32":  # Windows lacks the resource module
        assert mem["max_rss_gb"] > 0


@pytest.mark.parametrize(
    "changed",
    [{"hardware": "NVIDIA A100"}, {"versions": {"ase": "2"}}, {"dtype": None}],
)
def test_check_resume_labels(changed: dict[str, object]) -> None:
    """Resumed work must share every run label of the session before it."""
    runtime = {"hardware": "NVIDIA H200", "dtype": "float32", "versions": {"ase": "1"}}
    hpc.check_resume_labels(runtime, runtime | {"load_time_sec": 0.0}, "ctx")
    with pytest.raises(ValueError, match=f"ctx: cannot resume .*{next(iter(changed))}"):
        hpc.check_resume_labels(runtime, runtime | changed, "ctx")


@pytest.mark.parametrize(
    ("run_metadata", "yaml_cost", "n_problems"),
    [
        (RUNNER_RUN, RUNNER_COST, 0),
        (RUNNER_RUN | {"run_time_sec": 5.04}, RUNNER_COST, 1),
        (RUNNER_RUN, {}, 1),  # publishable artifact cost the YAML never registered
        # mixed-hardware merge kept tag/RSS only: nothing publishable, nothing missed
        ({"cost_source": "runner", "max_rss_gb": 1.0}, {}, 0),
        ({}, RUNNER_COST, 1),
        ({}, AUTHOR_COST, 0),  # YAML-only author report
        # author-run artifacts must match too, but their curated tag is not compared
        (RUNNER_RUN, RUNNER_COST | {"source": "author_reported"}, 0),
        (RUNNER_RUN, AUTHOR_COST, 2),  # hardware + run_time_sec differ
        # legacy artifacts are compared on their trustworthy keys only
        ({"hardware": "NVIDIA H200"}, RUNNER_COST | {"dtype": "float32"}, 1),
        (
            {"hardware": "NVIDIA H200", "run_time_sec": 5.0},
            RUNNER_COST | {"dtype": "float32"},
            0,
        ),
    ],
)
def test_cost_provenance_problems(
    tmp_path: Path,
    run_metadata: dict[str, object],
    yaml_cost: dict[str, object],
    n_problems: int,
) -> None:
    """YAML cost mappings must equal their artifact's recorded cost."""
    rel_path = "models/x/y/2026-01-01-discovery-run-info.json"
    (tmp_path / rel_path).parent.mkdir(parents=True)
    (tmp_path / rel_path).write_text(json.dumps({"run_metadata": run_metadata}))
    section = {"run_info_file": {"name": rel_path}} | (
        {"cost": yaml_cost} if yaml_cost else {}
    )
    problems = hpc.cost_provenance_problems({"discovery": section}, str(tmp_path))
    assert len(problems) == n_problems, problems
