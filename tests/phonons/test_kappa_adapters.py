"""Focused behavior tests for custom unified kappa adapters."""

from __future__ import annotations

import sys
from types import ModuleType, SimpleNamespace
from typing import TYPE_CHECKING, Any, cast
from unittest.mock import Mock

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator
from phonopy.structure.atoms import PhonopyAtoms

from matbench_discovery.enums import Key
from matbench_discovery.phonons.adapters.equflash import EquFlashKappaAdapter
from matbench_discovery.phonons.adapters.pet import pet_batch_forces
from matbench_discovery.phonons.adapters.standard import (
    FairchemKappaAdapter,
    StandardKappaAdapter,
)
from matbench_discovery.phonons.pipeline import KappaSettings
from matbench_discovery.phonons.thermal_conductivity import batched_displacement_forces

if TYPE_CHECKING:
    from collections.abc import Sequence

    from phono3py.api_phono3py import Phono3py


def displaced_h2(position_value: float) -> PhonopyAtoms:
    """Two-atom displaced supercell with every position coordinate = position_value."""
    return PhonopyAtoms(
        symbols=["H", "H"], cell=np.eye(3), positions=np.full((2, 3), position_value)
    )


class BatchCalculatorStub:
    """Echo atom positions as forces and record PET batch sizes."""

    def __init__(self) -> None:
        """Initialize an empty batch-size log."""
        self.batch_sizes: list[int] = []

    def compute_energy(
        self, atoms: list[Atoms], *, compute_forces_and_stresses: bool
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


class TensorStub(np.ndarray):
    """NumPy stand-in for the torch-geometric batch and torch tensor APIs EquFlash
    uses: a batch of position graphs is their concatenation.
    """

    @classmethod
    def from_data_list(cls, graphs: list[np.ndarray]) -> TensorStub:
        """Concatenate one chunk of graph stand-ins into a batch."""
        return np.concatenate(graphs).view(cls)

    def to(self, device: str) -> TensorStub:
        """Require the CPU device the test environment resolves to."""
        assert device == "cpu"
        return self

    def detach(self) -> TensorStub:
        """Return this test tensor."""
        return self

    def cpu(self) -> TensorStub:
        """Return this already-CPU test tensor."""
        return self

    def numpy(self) -> np.ndarray:
        """Return the underlying NumPy array."""
        return np.asarray(self)


class EquFlashTrainerStub:
    """Echo graph positions as forces while recording prediction batch sizes."""

    device: str | None = None

    def __init__(self) -> None:
        """Initialize an empty prediction log."""
        self.batch_sizes: list[int] = []

    def predict(
        self, graph_batch: TensorStub, *, per_image: bool, disable_tqdm: bool
    ) -> dict[str, TensorStub]:
        """Return the batch's concatenated graph positions as deterministic forces."""
        assert per_image is False
        assert disable_tqdm is True
        self.batch_sizes.append(len(graph_batch) // 2)  # 2-atom supercells
        return {"forces": graph_batch}


class EchoCalculator(Calculator):
    """ASE-only calculator (no batch method) returning zero energy, atom positions as
    forces and a fixed Voigt stress.
    """

    implemented_properties = ["energy", "forces", "stress"]  # noqa: RUF012

    def __init__(self, stress: Sequence[float] = (0,) * 6) -> None:
        """Store the Voigt stress (default zero) returned by every calculation."""
        super().__init__()
        self.stress = np.asarray(stress, dtype=float)

    def calculate(
        self,
        atoms: Atoms | None = None,
        properties: Sequence[str] = ("energy", "forces", "stress"),
        system_changes: Sequence[str] = (),
    ) -> None:
        """Populate deterministic results for the current structure."""
        super().calculate(atoms, properties, system_changes)
        if atoms is None:
            raise ValueError("atoms are required")
        self.results = {
            "energy": 0.0,
            "forces": atoms.positions.copy(),
            "stress": self.stress,
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

    def __init__(self, atoms: Atoms, *, symprec: float) -> None:
        """Retain constrained atoms and verify the relaxation tolerance."""
        self.atoms = atoms
        assert symprec == 0.01


def patch_relaxation(monkeypatch: pytest.MonkeyPatch, filter_class: object) -> None:
    """Install the cell filter (class or None), a no-op optimizer, constraint stub and
    initial, constrained and unconstrained space groups 225, 225, 221.
    """
    from matbench_discovery.phonons.adapters import standard

    monkeypatch.setattr(standard, "resolve_cell_filter", lambda _name: filter_class)
    monkeypatch.setattr(standard, "resolve_optimizer", lambda _name: NoOpOptimizer)
    monkeypatch.setattr(
        standard, "space_group_number", Mock(side_effect=[225, 225, 221])
    )
    monkeypatch.setattr(standard, "FixSymmetry", ConstraintStub)


@pytest.mark.parametrize(
    ("calculator_cls", "batch_size", "max_evaluations", "expected_values"),
    [
        (BatchCalculatorStub, 1, 1, (1, 0, 0)),
        (BatchCalculatorStub, 2, None, (1, 0, 3)),
        # current metatomic SymmetrizedCalculator exposes only the ASE API
        (EchoCalculator, 1, 1, (1, 0, 0)),
        (EchoCalculator, 2, None, (1, 0, 3)),
    ],
    ids=["batch_capped", "batch", "ase_capped", "ase"],
)
def test_pet_batched_force_set(
    calculator_cls: type[BatchCalculatorStub | EchoCalculator],
    batch_size: int,
    max_evaluations: int | None,
    expected_values: tuple[float, float, float],
) -> None:
    """PET evaluates up to the cap through its batch or ASE-only API and zero-fills
    null and deferred displacements.
    """
    calculator = calculator_cls()
    force_set = batched_displacement_forces(
        [displaced_h2(1), None, displaced_h2(3)],
        pet_batch_forces(cast("Calculator", calculator)),
        batch_size=batch_size,
        n_atoms=2,
        max_evaluations=max_evaluations,
    )
    expected = np.repeat(expected_values, 6).reshape(3, 2, 3).astype(float)
    np.testing.assert_array_equal(force_set, expected)
    if isinstance(calculator, BatchCalculatorStub):
        assert calculator.batch_sizes == [1 if max_evaluations else 2]


@pytest.mark.parametrize(
    ("enforce_relax_symm", "expected_filter_kwargs"),
    [
        (True, {"mask": [True, True, True, False, False, False]}),
        (False, {"mask": None}),
    ],
)
def test_standard_relax_only_masks_tilts_with_symmetry_constraint(
    monkeypatch: pytest.MonkeyPatch,
    enforce_relax_symm: bool,
    expected_filter_kwargs: dict[str, object],
) -> None:
    """Unconstrained legacy relaxations keep all cell degrees of freedom."""
    filter_calls: list[dict[str, object]] = []

    def filter_stub(atoms: Atoms, **kwargs: object) -> Atoms:
        """Capture filter options while exposing the original atoms."""
        filter_calls.append(kwargs)
        return atoms

    patch_relaxation(monkeypatch, filter_stub)
    atoms = Atoms("H", cell=np.eye(3), positions=[[0, 0, 0]], pbc=True)

    result = StandardKappaAdapter().relax(
        atoms,
        # compressive stress must report its magnitude, not a signed max of ~0
        EchoCalculator(stress=(-3, -3, -3, 0, 0, 0)),
        KappaSettings(enforce_relax_symm=enforce_relax_symm, max_steps=1),
        log_file=None,
    )

    assert filter_calls == [expected_filter_kwargs]
    assert result.max_stress is not None
    np.testing.assert_allclose(result.max_stress, [3, 0], rtol=0, atol=1e-12)


@pytest.mark.parametrize("enforce_relax_symm", [True, False])
def test_standard_relax_with_emt(enforce_relax_symm: bool) -> None:
    """Strained, rattled Cu relaxes under EMT; the symmetry constraint keeps Fm-3m,
    and relaxation_mode="none" returns the input untouched.
    """
    from ase.build import bulk
    from ase.calculators.emt import EMT

    atoms = bulk("Cu", "fcc", a=3.61 * 1.02, cubic=True)
    if not enforce_relax_symm:
        atoms.rattle(stdev=0.01, seed=0)
    settings = KappaSettings(
        enforce_relax_symm=enforce_relax_symm, force_max=1e-3, max_steps=300
    )
    result = StandardKappaAdapter().relax(atoms, EMT(), settings, log_file=None)
    if enforce_relax_symm:
        assert (result.initial_spg_num, result.final_spg_num) == (225, 225)
    assert 0 < result.n_steps < 300
    assert not result.reached_max_steps
    assert result.max_stress is not None
    assert np.abs(result.max_stress).max() < 1e-3
    assert len(result.atoms) == 4
    assert result.atoms.calc is None
    assert not result.atoms.constraints

    unrelaxed = bulk("Cu", "fcc", a=3.61 * 1.02, cubic=True)
    skipped = StandardKappaAdapter().relax(
        unrelaxed, EMT(), KappaSettings(relaxation_mode="none"), log_file=None
    )
    assert skipped.atoms is unrelaxed
    assert (skipped.n_steps, skipped.max_stress) == (0, None)


def test_fairchem_adapter_disables_training_scaler() -> None:
    """Fairchem preparation removes checkpoint force scaling exactly once."""
    calculator = SimpleNamespace(trainer=SimpleNamespace(scaler=object()))
    prepared = FairchemKappaAdapter().prepare_calculator(
        cast("Calculator", calculator), KappaSettings()
    )
    assert prepared is calculator
    assert calculator.trainer.scaler is None


def test_init_phono3py_forwards_q_point_mesh_and_plusminus(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Init forwards q_point_mesh and the displacement policy, rejects legacy q_mesh."""
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
    with pytest.raises(KeyError, match="q_point_mesh"):
        StandardKappaAdapter().init_phono3py(atoms, KappaSettings())
    atoms.info["q_point_mesh"] = atoms.info.pop("q_mesh")
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
    data_module.__dict__["Batch"] = TensorStub
    torch_module = ModuleType("torch")
    torch_module.__dict__["cuda"] = SimpleNamespace(is_available=bool)
    for module_name, module in (
        ("fairchem", ModuleType("fairchem")),
        ("fairchem.core", ModuleType("fairchem.core")),
        ("fairchem.core.preprocessing", preprocessing_module),
        ("torch", torch_module),
        ("torch_geometric", ModuleType("torch_geometric")),
        ("torch_geometric.data", data_module),
    ):
        monkeypatch.setitem(sys.modules, module_name, module)

    trainer = EquFlashTrainerStub()
    calculator = cast("Calculator", SimpleNamespace(trainer=trainer))
    phono3py = cast(
        "Phono3py",
        SimpleNamespace(
            supercell=Atoms("HH"),
            supercells_with_displacements=[
                displaced_h2(1),
                None,
                displaced_h2(3),
                displaced_h2(4),
            ],
            forces=None,
        ),
    )
    force_set = EquFlashKappaAdapter().calculate_fc3(
        phono3py, calculator, KappaSettings(batch_size=2), max_evaluations=2
    )

    assert trainer.batch_sizes == [2]
    expected = np.repeat([1, 0, 3, 0], 6).reshape(4, 2, 3).astype(float)
    np.testing.assert_array_equal(force_set, expected)
    assert phono3py.forces is force_set


def test_equflash_two_stage_relax_redirects_symmetry(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """EquFlash retains constrained geometry when unconstrained symmetry changes."""
    patch_relaxation(monkeypatch, None)
    atoms = Atoms("HH", cell=np.eye(3), positions=[[0] * 3, [0.5] * 3], pbc=True)
    result = EquFlashKappaAdapter().relax(
        atoms,
        EchoCalculator(),
        KappaSettings(relaxation_mode="two-stage"),
        log_file=None,
    )
    assert result.initial_spg_num == 225
    assert result.final_spg_num == 225
    assert result.redirected_to_symmetry is True
    assert result.n_steps == 2
