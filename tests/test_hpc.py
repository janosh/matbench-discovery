"""Tests for high-performance computing utilities."""

import os
from collections.abc import Callable, Sequence
from unittest.mock import mock_open, patch

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk
from pymatgen.core import Lattice, Structure

from matbench_discovery import hpc
from matbench_discovery.data import ase_atoms_from_zip
from matbench_discovery.enums import DataFiles


def make_ase_atoms(n_atoms: int) -> Atoms:
    """Helper to create ASE Atoms with n atoms."""
    return bulk("Cu") * (n_atoms, 1, 1)


def make_pmg_structure(n_atoms: int) -> Structure:
    """Helper to create pymatgen Structure with n atoms."""
    return Structure(
        Lattice.cubic(3.6),
        ["Cu"] * n_atoms,
        [[0, 0, 0]] * n_atoms,
        coords_are_cartesian=False,
        site_properties={"index": [0] * n_atoms},  # initialize site_properties
    )


@pytest.mark.parametrize("make_obj", [make_ase_atoms, make_pmg_structure])
def test_varying_size_objects(make_obj: Callable[[int], Atoms | Structure]) -> None:
    """Test chunking with structures of varying sizes."""
    # Create objects with 1, 8, 27, and 64 atoms
    objects = [make_obj(i) for i in range(1, 5)]
    chunks = hpc.chunk_by_lens(objects, n_chunks=2)

    # Calculate total atoms in each chunk
    atoms_per_chunk = [sum(len(obj) for obj in chunk) for chunk in chunks]

    # Largest structures should be split first
    assert len(chunks) == 2
    # Check if the difference in total atoms between chunks is minimal
    assert abs(atoms_per_chunk[0] - atoms_per_chunk[1]) <= max(
        len(obj) for obj in objects
    )


def test_empty_input() -> None:
    """Test chunking with empty input list."""
    assert hpc.chunk_by_lens([], n_chunks=1) == []


def test_invalid_n_chunks() -> None:
    """Test chunking with invalid number of chunks."""
    structures = [bulk("Cu")] * 3
    with pytest.raises(
        ValueError, match="n_chunks or chunk_size must be positive integer"
    ):
        hpc.chunk_by_lens(structures, n_chunks=0)


def test_n_chunks_larger_than_inputs() -> None:
    """Test when n_chunks is larger than number of inputs."""
    structures = [bulk("Cu")] * 3
    chunks = hpc.chunk_by_lens(structures, n_chunks=5)
    assert len(chunks) == len(structures)  # should reduce n_chunks to len(inputs)
    assert all(
        len(chunk) == 1 for chunk in chunks
    )  # each chunk should have 1 structure


def test_missing_chunk_params() -> None:
    """Test that error is raised when neither n_chunks nor chunk_size is specified."""
    structures = [bulk("Cu")] * 3
    with pytest.raises(
        ValueError, match="n_chunks or chunk_size must be positive integer"
    ):
        hpc.chunk_by_lens(structures)


def test_both_chunk_params() -> None:
    """Test that error is raised when both n_chunks and chunk_size are specified."""
    structures = [bulk("Cu")] * 3
    with pytest.raises(ValueError, match="Cannot specify both n_chunks and chunk_size"):
        hpc.chunk_by_lens(structures, n_chunks=2, chunk_size=10)


@pytest.mark.parametrize("make_obj", [make_ase_atoms, make_pmg_structure])
def test_chunk_size_basic(make_obj: Callable[[int], Atoms | Structure]) -> None:
    """Test basic chunk_size functionality."""
    # Create objects with 1, 8, 27, and 64 atoms
    objects = [make_obj(i) for i in range(1, 5)]
    total_atoms = sum(len(obj) for obj in objects)
    chunk_size = total_atoms // 3  # should create ~3 chunks

    chunks = hpc.chunk_by_lens(objects, chunk_size=chunk_size)
    atoms_per_chunk = [sum(len(obj) for obj in chunk) for chunk in chunks]

    # Check that most chunks are close to the target size
    assert all(size <= chunk_size * 1.5 for size in atoms_per_chunk)
    assert sum(atoms_per_chunk) == total_atoms


@pytest.mark.parametrize("make_obj", [make_ase_atoms, make_pmg_structure])
def test_chunk_size_small(make_obj: Callable[[int], Atoms | Structure]) -> None:
    """Test chunking with chunk_size smaller than smallest object."""
    objects = [make_obj(5) for _ in range(3)]  # 3 objects with 5 atoms each
    chunks = hpc.chunk_by_lens(objects, chunk_size=3)

    # Should create one chunk per object since each object is larger than chunk_size
    assert len(chunks) == 5
    assert [len(chunk) for chunk in chunks] == [1, 1, 1, 0, 0]


@pytest.mark.parametrize("make_obj", [make_ase_atoms, make_pmg_structure])
def test_chunk_size_large(make_obj: Callable[[int], Atoms | Structure]) -> None:
    """Test chunking with chunk_size larger than total size."""
    objects = [make_obj(i) for i in range(1, 4)]  # objects with 1, 2, 3 atoms
    chunks = hpc.chunk_by_lens(objects, chunk_size=10)

    # Should create a single chunk containing all objects
    assert len(chunks) == 1
    assert len(chunks[0]) == 3


@pytest.mark.parametrize("make_obj", [make_ase_atoms, make_pmg_structure])
def test_equal_size_objects(make_obj: Callable[[int], Atoms | Structure]) -> None:
    """Test chunking with objects of equal size."""
    objects = [make_obj(1) for _ in range(6)]  # each has 1 atom
    chunks = hpc.chunk_by_lens(objects, n_chunks=2)
    assert len(chunks) == 2
    assert len(chunks[0]) == len(chunks[1]) == 3  # should split evenly
    assert sum(len(obj) for obj in chunks[0]) == sum(len(obj) for obj in chunks[1])


def test_different_sized_objects() -> None:
    """Test chunking with different types of sized objects."""
    # Mix of lists and tuples of different lengths
    inputs: Sequence[Sequence[int]] = [(1, 2, 3), [1], (1, 2), [1, 2, 3, 4]]
    chunks = hpc.chunk_by_lens(inputs, n_chunks=2)

    # Calculate total length in each chunk
    lens_per_chunk = [sum(len(obj) for obj in chunk) for chunk in chunks]

    assert len(chunks) == 2
    # Check if the difference in total lengths between chunks is minimal
    assert abs(lens_per_chunk[0] - lens_per_chunk[1]) <= max(len(obj) for obj in inputs)


@pytest.mark.parametrize("make_obj", [make_ase_atoms, make_pmg_structure])
def test_chunk_size_distribution(make_obj: Callable[[int], Atoms | Structure]) -> None:
    """Test the distribution of chunk sizes."""
    # Create objects with exponentially increasing sizes: 1, 8, 27, 64, 125 atoms
    objects = [make_obj(i) for i in range(1, 6)]
    chunks = hpc.chunk_by_lens(objects, n_chunks=3)

    # Calculate atoms per chunk
    atoms_per_chunk = np.array([sum(len(obj) for obj in chunk) for chunk in chunks])
    mean_atoms = atoms_per_chunk.mean()

    # Check that no chunk deviates too much from the mean
    max_deviation = abs(atoms_per_chunk - mean_atoms).max()
    total_atoms = sum(len(obj) for obj in objects)
    assert max_deviation < total_atoms / len(
        chunks
    )  # deviation should be less than average chunk size


@pytest.mark.parametrize("make_obj", [make_ase_atoms, make_pmg_structure])
def test_order_preservation(make_obj: Callable[[int], Atoms | Structure]) -> None:
    """Test that relative order of equal-sized structures is preserved within chunks."""
    # Create 6 identical objects but with different info
    objects = [make_obj(1) for _ in range(6)]
    for idx, obj in enumerate(objects):
        if hasattr(obj, "info"):
            obj.info["index"] = idx  # ASE Atoms
        else:
            obj.site_properties["index"] = [idx] * len(obj)  # pymatgen Structure

    chunks = hpc.chunk_by_lens(objects, n_chunks=2)

    # Within each chunk, indices should be monotonically decreasing
    for chunk in chunks:
        indices = []
        for obj in chunk:
            if hasattr(obj, "info"):
                indices.append(obj.info["index"])  # ASE Atoms
            else:
                indices.append(obj.site_properties["index"][0])  # pymatgen Structure
        assert indices == sorted(indices, reverse=True)


@pytest.mark.parametrize("make_obj", [make_ase_atoms, make_pmg_structure])
def test_reproducibility(make_obj: Callable[[int], Atoms | Structure]) -> None:
    """Test that chunking is deterministic."""
    objects = [make_obj(i) for i in range(1, 5)]
    chunks1 = hpc.chunk_by_lens(objects, n_chunks=2)
    chunks2 = hpc.chunk_by_lens(objects, n_chunks=2)

    # Compare object identities in each chunk
    for chunk1, chunk2 in zip(chunks1, chunks2, strict=True):
        assert len(chunk1) == len(chunk2)
        for obj1, obj2 in zip(chunk1, chunk2, strict=True):
            assert obj1 is obj2  # should be the exact same object


@pytest.mark.slow
def test_wbm_chunking() -> None:
    """Test chunking with WBM initial atoms dataset."""
    atoms_list = ase_atoms_from_zip(DataFiles.wbm_initial_atoms.path)
    n_chunks = 100

    # Verify distribution of atoms per structure
    atoms_per_struct = np.array([len(atoms) for atoms in atoms_list])
    unique_sizes, counts = np.unique(atoms_per_struct, return_counts=True)
    expected_sizes = (
        {2: 1596, 3: 6040, 4: 45937, 5: 18332, 6: 35209, 7: 8691, 8: 38918}
        | {9: 21670, 10: 32982, 11: 9097, 12: 36499, 15: 22, 16: 429, 20: 837, 24: 60}
        | {25: 1, 30: 164, 32: 126, 40: 195, 45: 2, 48: 91, 60: 31, 64: 11, 70: 5}
        | {80: 11, 90: 3, 100: 4}
    )
    np.testing.assert_array_equal(unique_sizes, list(expected_sizes))
    assert sorted(counts) == sorted(expected_sizes.values())
    assert len(atoms_list) == 256963, "WBM dataset should have 256,963 structures"

    # Split into chunks
    chunks = hpc.chunk_by_lens(atoms_list, n_chunks=n_chunks)

    # Verify chunk statistics
    atoms_per_chunk = np.array([sum(len(atoms) for atoms in chunk) for chunk in chunks])
    structs_per_chunk = np.array([len(chunk) for chunk in chunks])

    # Since all structures have same size, chunks should be nearly equal
    assert atoms_per_chunk.mean() == pytest.approx(19941.3, rel=1e-3)
    assert atoms_per_chunk.std() == pytest.approx(0.906, rel=1e-3)
    assert atoms_per_chunk.min() == 19940.0
    assert atoms_per_chunk.max() == 19942.0

    assert structs_per_chunk.mean() == pytest.approx(2569.6, rel=1e-3)
    assert structs_per_chunk.std() == pytest.approx(0.4828, rel=1e-3)
    assert structs_per_chunk.min() == 2569
    assert structs_per_chunk.max() == 2570

    # Verify relative standard deviations are very small
    atoms_rsd = 100 * atoms_per_chunk.std() / atoms_per_chunk.mean()
    structs_rsd = 100 * structs_per_chunk.std() / structs_per_chunk.mean()
    assert atoms_rsd < 0.1, "Relative std dev in atoms per chunk should be < 0.1%"
    assert structs_rsd < 0.1, (
        "Relative std dev in structures per chunk should be < 0.1%"
    )


@pytest.mark.parametrize("py_file_path", [None, "path/to/file.py"])
@pytest.mark.parametrize("partition", [None, "fake-partition"])
@pytest.mark.parametrize("time", [None, "0:0:1"])
@pytest.mark.parametrize("account", [None, "fake-account"])
@pytest.mark.parametrize("pre_cmd", [None, "module load pytorch;", "ENV_VAR=42"])
@pytest.mark.parametrize("submit_as_temp_file", [True, False])
def test_slurm_submit(
    capsys: pytest.CaptureFixture[str],
    py_file_path: str | None,
    partition: str | None,
    time: str | None,
    account: str | None,
    pre_cmd: str | None,
    submit_as_temp_file: bool,
) -> None:
    job_name = "test_job"
    out_dir = "tmp"

    kwargs = dict(
        job_name=job_name,
        out_dir=out_dir,
        time=time,
        partition=partition,
        account=account,
        py_file_path=py_file_path,
        slurm_flags="--foo",
        pre_cmd=pre_cmd,
        submit_as_temp_file=submit_as_temp_file,
    )

    hpc.slurm_submit(**kwargs)  # type: ignore[arg-type]

    stdout, stderr = capsys.readouterr()
    # check slurm_submit() did nothing in normal mode
    assert stdout == stderr == ""

    with patch.dict(os.environ, {"SLURM_JOB_ID": "1234"}, clear=True):
        slurm_vars = hpc.slurm_submit(**kwargs)  # type: ignore[arg-type]
    expected_slurm_vars = dict(slurm_job_id="1234", slurm_flags="--foo")
    if time is not None:
        expected_slurm_vars["slurm_timelimit"] = time
    if pre_cmd and not pre_cmd.strip().endswith(";"):
        pre_cmd += ";"
    if pre_cmd is not None and pre_cmd != "":
        expected_slurm_vars["pre_cmd"] = pre_cmd
    assert slurm_vars == expected_slurm_vars

    # check slurm_submit() prints cmd and calls subprocess.run() in submit mode
    with (
        pytest.raises(SystemExit),
        patch("sys.argv", ["slurm-submit"]),
        patch("matbench_discovery.hpc.subprocess.run") as mock_subprocess_run,
        patch(
            "matbench_discovery.hpc.tempfile.mkdtemp",
            return_value="/tmp/slurm_job_123",
        ),
        patch("matbench_discovery.hpc.shutil.copy2") as mock_copy2,
        patch("builtins.open", mock_open()),
        patch.dict(os.environ, {"SLURM_JOB_ID": "1234"}, clear=True),
    ):
        hpc.slurm_submit(**kwargs)  # type: ignore[arg-type]

    assert mock_subprocess_run.call_count == 1

    expected_py_file_path = py_file_path or __file__
    if submit_as_temp_file:
        expected_py_file_path = (
            f"/tmp/slurm_job_123/{os.path.basename(expected_py_file_path)}"
        )
        assert mock_copy2.called
        assert mock_copy2.call_args[0][0] == (py_file_path or __file__)
        assert mock_copy2.call_args[0][1] == expected_py_file_path
    else:
        assert not mock_copy2.called

    sbatch_cmd = (
        f"sbatch --job-name {job_name} --output {out_dir}/slurm-%A.log --foo "
        f"--wrap {pre_cmd + ' ' if pre_cmd else ''}python {expected_py_file_path}"
    ).replace(" --", "\n  --")
    for flag in (f"{time=!s}", f"{account=!s}", f"{partition=!s}"):
        key, val = flag.split("=")
        if val != "None":
            sbatch_cmd += f"\n  --{key} {val}"

    stdout, stderr = capsys.readouterr()
    assert sbatch_cmd in stdout
    assert stderr == ""


def test_get_calling_file_path() -> None:
    assert hpc._get_calling_file_path(frame=1) == __file__  # noqa: SLF001

    def wrapper(frame: int) -> str:
        return hpc._get_calling_file_path(frame)  # noqa: SLF001

    assert wrapper(frame=2) == __file__
