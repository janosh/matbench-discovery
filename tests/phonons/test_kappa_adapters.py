"""Focused behavior tests for custom unified kappa adapters."""

from __future__ import annotations

import sys
from types import ModuleType, SimpleNamespace
from typing import TYPE_CHECKING, Any, cast

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator
from pymatviz.enums import Key

from matbench_discovery.phonons.adapters.equflash import EquFlashKappaAdapter
from matbench_discovery.phonons.adapters.pet import _batched_force_set
from matbench_discovery.phonons.adapters.standard import (
    FairchemKappaAdapter,
    StandardKappaAdapter,
)
from matbench_discovery.phonons.pipeline import KappaSettings

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

    from phono3py.api_phono3py import Phono3py
    from phonopy.structure.atoms import PhonopyAtoms


class FakePhonopyAtoms:
    """Minimal displaced-supercell shape consumed by PET batching."""

    def __init__(self, position_value: float) -> None:
        """Create a two-atom periodic supercell with constant positions."""
        self.symbols = ["H", "H"]
        self.cell = np.eye(3)
        self.positions = np.full((2, 3), position_value)


class BatchCalculatorStub:
    """Echo atom positions as forces and record PET batch sizes."""

    def __init__(self) -> None:
        """Initialize an empty batch-size log."""
        self.batch_sizes: list[int] = []

    def compute_energy(
        self,
        atoms: list[Atoms],
        *,
        compute_forces_and_stresses: bool,
    ) -> dict[str, Any]:
        """Return positions as deterministic forces for one requested batch."""
        assert compute_forces_and_stresses is True
        self.batch_sizes.append(len(atoms))
        return {"forces": np.asarray([item.positions for item in atoms])}


class GraphConverterStub:
    """Represent each EquFlash graph by its source atom positions."""

    def __init__(self, *, r_edges: bool) -> None:
        """Require the same edge option as the production converter."""
        assert r_edges is False

    def convert(self, atoms: Atoms) -> np.ndarray:
        """Return positions as a lightweight graph stand-in."""
        return atoms.positions.copy()


class GraphBatchStub:
    """Minimal torch-geometric batch carrying converted position arrays."""

    def __init__(self, graphs: list[np.ndarray]) -> None:
        """Store one production-sized graph chunk."""
        self.graphs = graphs

    @classmethod
    def from_data_list(cls, graphs: list[np.ndarray]) -> GraphBatchStub:
        """Build a batch from one list of graph stand-ins."""
        return cls(graphs)

    def to(self, device: str) -> GraphBatchStub:
        """Record no device state and return this CPU test batch."""
        assert device == "cpu"
        return self


class ArrayTensorStub:
    """Expose the tensor conversion chain used by the EquFlash adapter."""

    def __init__(self, array: np.ndarray) -> None:
        """Store one force array."""
        self.array = array

    def detach(self) -> ArrayTensorStub:
        """Return this immutable test tensor."""
        return self

    def cpu(self) -> ArrayTensorStub:
        """Return this already-CPU test tensor."""
        return self

    def numpy(self) -> np.ndarray:
        """Return the underlying NumPy force array."""
        return self.array


class EquFlashTrainerStub:
    """Echo graph positions as forces while recording prediction batch sizes."""

    device: str | None = None

    def __init__(self) -> None:
        """Initialize an empty prediction log."""
        self.batch_sizes: list[int] = []

    def predict(
        self,
        graph_batch: GraphBatchStub,
        *,
        per_image: bool,
        disable_tqdm: bool,
    ) -> dict[str, ArrayTensorStub]:
        """Return concatenated graph positions as deterministic forces."""
        assert per_image is False
        assert disable_tqdm is True
        self.batch_sizes.append(len(graph_batch.graphs))
        return {"forces": ArrayTensorStub(np.concatenate(graph_batch.graphs, axis=0))}


class ZeroStressCalculator(Calculator):
    """ASE calculator returning zero energy, forces, and stress."""

    implemented_properties = ["energy", "forces", "stress"]  # noqa: RUF012

    def calculate(
        self,
        atoms: Atoms | None = None,
        properties: Sequence[str] = ("energy", "forces", "stress"),
        system_changes: Sequence[str] = (),
    ) -> None:
        """Populate zero-valued results for relaxation tests."""
        super().calculate(atoms, properties, system_changes)
        if atoms is None:
            raise ValueError("atoms are required")
        self.results = {
            "energy": 0.0,
            "forces": np.zeros((len(atoms), 3)),
            "stress": np.zeros(6),
        }


class NoOpOptimizer:
    """ASE optimizer stub that reports one completed step."""

    def __init__(self, atoms: object, logfile: str | None = None) -> None:
        """Store constructor arguments without changing the structure."""
        self.atoms = atoms
        self.logfile = logfile
        self.n_steps = 0

    def run(self, fmax: float, steps: int) -> None:
        """Record one step for any non-empty optimization budget."""
        del fmax
        self.n_steps = min(steps, 1)

    def get_number_of_steps(self) -> int:
        """Return the recorded optimizer step count."""
        return self.n_steps


class ConstraintStub:
    """Copyable placeholder for ASE's symmetry constraint."""

    def __init__(self, atoms: Atoms) -> None:
        """Retain the constrained atoms only for debugging."""
        self.atoms = atoms


class SpaceGroupSequence:
    """Return deterministic initial, constrained, and unconstrained groups."""

    def __init__(self) -> None:
        """Initialize the expected three-call sequence."""
        self.values = iter((225, 225, 221))

    def __call__(self, atoms: Atoms, symprec: float) -> int:
        """Return the next test space group."""
        del atoms, symprec
        return next(self.values)


def resolve_no_filter(name: str | None) -> None:
    """Return no ASE filter for the two-stage relaxation unit test."""
    del name


def resolve_noop_optimizer(name: str) -> type[NoOpOptimizer]:
    """Return the no-op optimizer used by the two-stage relaxation test."""
    del name
    return NoOpOptimizer


def patch_relaxation(
    monkeypatch: pytest.MonkeyPatch,
    filter_resolver: Callable[[str | None], object],
) -> None:
    """Install deterministic optimizer, symmetry, filter, and constraint stubs."""
    from matbench_discovery.phonons.adapters import standard

    monkeypatch.setattr(standard, "resolve_cell_filter", filter_resolver)
    monkeypatch.setattr(standard, "resolve_optimizer", resolve_noop_optimizer)
    monkeypatch.setattr(standard, "space_group_number", SpaceGroupSequence())
    monkeypatch.setattr(standard, "FixSymmetry", ConstraintStub)


def test_pet_batched_force_set_preserves_shape_when_capped() -> None:
    """PET evaluates the cap and zero-fills null and deferred displacements."""
    calculator = BatchCalculatorStub()
    displacements = cast(
        "list[PhonopyAtoms | None]",
        [FakePhonopyAtoms(1), None, FakePhonopyAtoms(3)],
    )
    force_set = _batched_force_set(
        displacements,
        cast("Calculator", calculator),
        batch_size=1,
        n_atoms=2,
        max_evaluations=1,
    )
    assert calculator.batch_sizes == [1]
    assert np.all(force_set[0] == 1)
    assert np.all(force_set[1] == 0)
    assert np.all(force_set[2] == 0)


def test_pet_force_set_supports_ase_only_symmetrizer_api() -> None:
    """PET falls back to ASE calls when its symmetrizer has no batch method."""
    displacements = cast(
        "list[PhonopyAtoms | None]",
        [FakePhonopyAtoms(1), None, FakePhonopyAtoms(3)],
    )
    force_set = _batched_force_set(
        displacements, ZeroStressCalculator(), batch_size=2, n_atoms=2
    )
    assert force_set.shape == (3, 2, 3)
    assert np.all(force_set == 0)


@pytest.mark.parametrize(
    ("enforce_relax_symm", "expected_filter_kwargs"),
    [(True, {"mask": [True, True, True, False, False, False]}), (False, {})],
)
def test_standard_relax_only_masks_tilts_with_symmetry_constraint(
    monkeypatch: pytest.MonkeyPatch,
    enforce_relax_symm: bool,
    expected_filter_kwargs: dict[str, object],
) -> None:
    """Unconstrained legacy relaxations keep all cell degrees of freedom."""
    from matbench_discovery.phonons.adapters import standard

    filter_calls: list[dict[str, object]] = []

    def filter_stub(atoms: Atoms, **kwargs: object) -> Atoms:
        """Capture filter options while exposing the original atoms."""
        filter_calls.append(kwargs)
        return atoms

    def resolve_filter_stub(name: str | None) -> Callable[..., Atoms]:
        """Return the capturing filter regardless of configured name."""
        del name
        return filter_stub

    patch_relaxation(monkeypatch, resolve_filter_stub)
    atoms = Atoms("H", cell=np.eye(3), positions=[[0, 0, 0]], pbc=True)

    standard.StandardKappaAdapter().relax(
        atoms,
        ZeroStressCalculator(),
        KappaSettings(enforce_relax_symm=enforce_relax_symm, max_steps=1),
        log_file=None,
    )

    assert filter_calls == [expected_filter_kwargs]


def test_fairchem_adapter_disables_training_scaler() -> None:
    """Fairchem preparation removes checkpoint force scaling exactly once."""
    calculator = SimpleNamespace(trainer=SimpleNamespace(scaler=object()))
    prepared = FairchemKappaAdapter().prepare_calculator(
        cast("Calculator", calculator), KappaSettings()
    )
    assert prepared is calculator
    assert calculator.trainer.scaler is None


def test_orb_adapter_accepts_q_mesh_and_plusminus(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """ORB forwards legacy q_mesh metadata and its displacement policy."""
    from matbench_discovery.phonons import thermal_conductivity

    initializer_kwargs: dict[str, object] = {}

    def init_phono3py_stub(atoms: Atoms, **kwargs: object) -> Phono3py:
        """Capture initialization arguments and return an opaque phono3py stub."""
        assert atoms.info[Key.mat_id] == "mp-orb"
        initializer_kwargs.update(kwargs)
        return cast("Phono3py", object())

    monkeypatch.setattr(thermal_conductivity, "init_phono3py", init_phono3py_stub)
    atoms = Atoms("H", cell=np.eye(3), positions=[[0, 0, 0]], pbc=True)
    atoms.info |= {
        Key.mat_id: "mp-orb",
        "fc2_supercell": np.eye(3, dtype=int),
        "fc3_supercell": np.eye(3, dtype=int),
        "q_mesh": [4, 5, 6],
    }
    StandardKappaAdapter().init_phono3py(atoms, KappaSettings(is_plusminus=True))
    assert initializer_kwargs["q_point_mesh"] == (4, 5, 6)
    assert initializer_kwargs["is_plusminus"] is True


def test_equflash_fc3_streams_graph_batches_and_preserves_null_entries(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """EquFlash converts at most one prediction batch and zero-fills null cells."""
    preprocessing_module = ModuleType("fairchem.core.preprocessing")
    preprocessing_module.__dict__["AtomsToGraphs"] = GraphConverterStub
    data_module = ModuleType("torch_geometric.data")
    data_module.__dict__["Batch"] = GraphBatchStub
    for module_name, module in (
        ("fairchem", ModuleType("fairchem")),
        ("fairchem.core", ModuleType("fairchem.core")),
        ("fairchem.core.preprocessing", preprocessing_module),
        ("torch_geometric", ModuleType("torch_geometric")),
        ("torch_geometric.data", data_module),
    ):
        monkeypatch.setitem(sys.modules, module_name, module)

    trainer = EquFlashTrainerStub()
    monkeypatch.setattr("torch.cuda.is_available", bool)
    calculator = cast("Calculator", SimpleNamespace(trainer=trainer))
    phono3py = cast(
        "Phono3py",
        SimpleNamespace(
            supercell=Atoms("HH"),
            supercells_with_displacements=[
                FakePhonopyAtoms(1),
                None,
                FakePhonopyAtoms(3),
                FakePhonopyAtoms(4),
            ],
            forces=None,
        ),
    )
    force_set = EquFlashKappaAdapter().calculate_fc3(
        phono3py,
        calculator,
        KappaSettings(batch_size=2),
        max_evaluations=2,
    )

    assert trainer.batch_sizes == [2]
    assert np.all(force_set[0] == 1)
    assert np.all(force_set[1] == 0)
    assert np.all(force_set[2] == 3)
    assert np.all(force_set[3] == 0)
    assert phono3py.forces is force_set


def test_equflash_two_stage_relax_redirects_symmetry(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """EquFlash retains constrained geometry when unconstrained symmetry changes."""
    patch_relaxation(monkeypatch, resolve_no_filter)
    atoms = Atoms(
        "HH",
        cell=np.eye(3),
        positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
        pbc=True,
    )
    result = EquFlashKappaAdapter().relax(
        atoms,
        ZeroStressCalculator(),
        KappaSettings(relaxation_mode="two-stage"),
        log_file=None,
    )
    assert result.initial_spg_num == 225
    assert result.final_spg_num == 225
    assert result.redirected_to_symmetry is True
    assert result.n_steps == 2
