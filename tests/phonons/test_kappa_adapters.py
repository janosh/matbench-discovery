"""Tests for the standard and backend-specific kappa adapters."""

from __future__ import annotations

import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import Any

import numpy as np
import pytest
from matscikit import EMT, Structure

from matbench_discovery.phonons import thermal_conductivity
from matbench_discovery.phonons.adapters import get_kappa_adapter
from matbench_discovery.phonons.adapters.equflash import (
    EquFlashBatchedModel,
    EquFlashKappaAdapter,
)
from matbench_discovery.phonons.adapters.pet import PetBatchedModel, PetKappaAdapter
from matbench_discovery.phonons.adapters.standard import (
    FairchemKappaAdapter,
    StandardKappaAdapter,
)
from matbench_discovery.phonons.pipeline import KappaSettings, RelaxationMode


def cu_cell(strain: float = 0.0, rattle: float = 0.0) -> Structure:
    """Conventional fcc Cu (4 atoms), optionally strained and rattled."""
    a = 3.61 * (1 + strain)
    frac = np.array([[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]], float)
    frac += rattle * np.random.default_rng(0).standard_normal(frac.shape)
    return Structure(a * np.eye(3), ["Cu"] * 4, frac, fractional=True)


@pytest.mark.parametrize(
    ("model_key", "adapter_cls"),
    [
        ("pet_oam_xl_1_0_0", PetKappaAdapter),
        ("equflash_29m_oam", EquFlashKappaAdapter),
        ("eqv2_m_oam_mp_salex", FairchemKappaAdapter),
        ("mace_mp_0", StandardKappaAdapter),
    ],
)
def test_get_kappa_adapter(model_key: str, adapter_cls: type) -> None:
    """Model keys dispatch to their backend-specific adapter."""
    assert type(get_kappa_adapter(model_key)) is adapter_cls


def test_fairchem_adapter_disables_training_scaler() -> None:
    """The wrapped calculator's trainer loses its training-only force scaler."""
    model = SimpleNamespace(
        calculator=SimpleNamespace(trainer=SimpleNamespace(scaler=object()))
    )
    prepared = FairchemKappaAdapter().prepare_model(model, KappaSettings())  # ty: ignore[invalid-argument-type]
    assert prepared is model
    assert model.calculator.trainer.scaler is None


def test_init_phonons_forwards_mesh_and_plusminus(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """PhononDB metadata and plus/minus settings reach the phonon setup."""
    calls: list[dict[str, object]] = []

    def init_stub(atoms: Structure, **kwargs: object) -> object:
        """Capture setup arguments."""
        del atoms
        calls.append(kwargs)
        return object()

    monkeypatch.setattr(thermal_conductivity, "init_phonons", init_stub)
    atoms = cu_cell()
    atoms.info |= {
        "fc2_supercell": np.diag([2, 2, 2]),
        "fc3_supercell": np.eye(3),
        "q_point_mesh": [3, 4, 5],
    }
    StandardKappaAdapter().init_phonons(atoms, KappaSettings())
    StandardKappaAdapter().init_phonons(atoms, KappaSettings(is_plusminus=True))
    assert [call["q_point_mesh"] for call in calls] == [(3, 4, 5)] * 2
    assert [call["is_plusminus"] for call in calls] == ["auto", True]


@pytest.mark.parametrize("enforce_relax_symm", [True, False])
def test_standard_relax_with_emt(enforce_relax_symm: bool) -> None:
    """Strained, rattled Cu relaxes; the symmetry constraint keeps its space group."""
    atoms = cu_cell(strain=0.02, rattle=0.0 if enforce_relax_symm else 0.01)
    settings = KappaSettings(
        enforce_relax_symm=enforce_relax_symm, force_max=1e-3, max_steps=300
    )
    result = StandardKappaAdapter().relax(atoms, EMT(), settings, log_file=None)
    assert result.initial_spg_num == (
        225 if enforce_relax_symm else result.initial_spg_num
    )
    if enforce_relax_symm:
        assert result.final_spg_num == 225
    assert 0 < result.n_steps < 300
    assert not result.reached_max_steps
    assert result.max_stress is not None
    assert np.abs(result.max_stress).max() < 1e-3
    assert len(result.atoms) == 4


def test_relax_skips_with_mode_none() -> None:
    """relaxation_mode='none' returns the input untouched without model calls."""
    atoms = cu_cell(strain=0.02)
    result = StandardKappaAdapter().relax(
        atoms, EMT(), KappaSettings(relaxation_mode="none"), log_file=None
    )
    assert result.atoms is atoms
    assert (result.n_steps, result.max_stress) == (0, None)


def test_two_stage_relax_writes_log(tmp_path: Path) -> None:
    """Two-stage relaxation keeps symmetry and logs both stages."""
    log_file = str(tmp_path / "logs" / "cu.log")
    result = EquFlashKappaAdapter().relax(
        cu_cell(strain=0.02),
        EMT(),
        KappaSettings(relaxation_mode="two-stage", force_max=1e-3),
        log_file=log_file,
    )
    assert result.final_spg_num == 225
    assert result.redirected_to_symmetry is False
    with open(log_file) as file:
        assert [line.split(":")[0] for line in file] == ["stage1", "stage2"]


def test_equflash_batched_model(monkeypatch: pytest.MonkeyPatch) -> None:
    """Forces-only batches run as one trainer call and split back per structure."""

    class GraphConverterStub:
        """Return each structure's positions as its 'graph'."""

        def __init__(self, *, r_edges: bool) -> None:
            assert r_edges is False

        def convert(self, atoms: Structure) -> np.ndarray:
            return atoms.positions

    class BatchStub:
        """Concatenate graphs; .to(device) is a no-op."""

        @staticmethod
        def from_data_list(graphs: list[np.ndarray]) -> SimpleNamespace:
            stacked = np.concatenate(graphs)
            return SimpleNamespace(to=lambda _device: stacked)

    predictions: list[int] = []

    def predict(batch: np.ndarray, **_kwargs: object) -> dict[str, Any]:
        predictions.append(len(batch))
        return {
            "forces": SimpleNamespace(
                detach=lambda: SimpleNamespace(
                    cpu=lambda: SimpleNamespace(numpy=lambda: -batch)
                )
            )
        }

    preprocessing = ModuleType("fairchem.core.preprocessing")
    preprocessing.__dict__["AtomsToGraphs"] = GraphConverterStub
    geometric = ModuleType("torch_geometric.data")
    geometric.__dict__["Batch"] = BatchStub
    monkeypatch.setitem(sys.modules, "fairchem.core.preprocessing", preprocessing)
    monkeypatch.setitem(sys.modules, "torch_geometric.data", geometric)
    pytest.importorskip("ase")  # structure_to_ase backs the graph conversion
    calculator = SimpleNamespace(trainer=SimpleNamespace(predict=predict, device="cpu"))
    model = EquFlashKappaAdapter().prepare_model(calculator, KappaSettings())  # ty: ignore[invalid-argument-type]
    assert isinstance(model, EquFlashBatchedModel)
    structures = [cu_cell(), cu_cell(strain=0.01)]
    results = model.compute(structures, ["forces"])
    assert predictions == [8]
    for structure, result in zip(structures, results, strict=True):
        np.testing.assert_allclose(result["forces"], -structure.positions)
    with pytest.raises(TypeError, match="no trainer"):
        EquFlashBatchedModel(SimpleNamespace())


def test_pet_batched_model_uses_compute_energy() -> None:
    """compute_energy batches energy+forces; other properties use the ASE path."""
    pytest.importorskip("ase")
    calls: list[int] = []

    def compute_energy(
        atoms_list: list[Any], *, compute_forces_and_stresses: bool
    ) -> dict:
        assert compute_forces_and_stresses
        calls.append(len(atoms_list))
        return {
            "energy": np.arange(len(atoms_list), dtype=float),
            "forces": np.array([atoms.positions for atoms in atoms_list]),
        }

    model = PetBatchedModel(SimpleNamespace(compute_energy=compute_energy))
    structures = [cu_cell(), cu_cell(strain=0.01), cu_cell(strain=0.02)]
    results = model.compute(structures, ["energy", "forces"])
    assert calls == [3]
    assert [result["energy"] for result in results] == [0.0, 1.0, 2.0]
    np.testing.assert_allclose(results[2]["forces"], structures[2].positions)


def spy_relax_structures(monkeypatch: pytest.MonkeyPatch) -> list[dict[str, Any]]:
    """Record relax_structures kwargs per stage and return the input unchanged."""
    from matbench_discovery.phonons.adapters import standard

    calls: list[dict[str, Any]] = []

    def relax_stub(
        structures: list[Structure], model: object, **kwargs: object
    ) -> list[dict[str, Any]]:
        """Capture options; report a converged zero-stress single step."""
        del model
        calls.append(kwargs)
        return [
            {"structure": s, "n_steps": 1, "stress": np.zeros((3, 3)), "error": None}
            for s in structures
        ]

    monkeypatch.setattr(standard, "relax_structures", relax_stub)
    return calls


@pytest.mark.parametrize(
    ("mode", "enforce_relax_symm", "expected"),
    [
        ("single-stage", True, [(0.01, (True, True, True, False, False, False))]),
        ("single-stage", False, [(None, None)]),
        (
            "two-stage",
            True,
            [
                (0.01, (True, True, True, False, False, False)),
                (None, (True, True, True, False, False, False)),
            ],
        ),
        ("two-stage", False, [(None, (True, True, True, False, False, False))] * 2),
    ],
)
def test_relax_stages_mask_tilts_and_fix_symmetry(
    monkeypatch: pytest.MonkeyPatch,
    mode: RelaxationMode,
    enforce_relax_symm: bool,
    expected: list[tuple],
) -> None:
    """No-tilt masks and the symmetry constraint reach exactly the intended stages."""
    calls = spy_relax_structures(monkeypatch)
    settings = KappaSettings(
        relaxation_mode=mode, enforce_relax_symm=enforce_relax_symm, relax_symprec=0.01
    )
    adapter = (
        StandardKappaAdapter() if mode == "single-stage" else EquFlashKappaAdapter()
    )
    adapter.relax(cu_cell(), EMT(), settings, log_file=None)
    assert [(call["fix_symmetry"], call["cell_mask"]) for call in calls] == expected
    assert all(call["stress"] is True for call in calls)


def test_two_stage_relax_redirects_to_constrained_symmetry(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """When the free stage changes the space group, the constrained geometry is kept."""
    from matbench_discovery.phonons.adapters import standard

    spy_relax_structures(monkeypatch)
    groups = iter([225, 225, 221])  # initial, stage 1, stage 2

    monkeypatch.setattr(standard, "space_group_number", lambda *_args: next(groups))
    result = EquFlashKappaAdapter().relax(
        cu_cell(), EMT(), KappaSettings(relaxation_mode="two-stage"), log_file=None
    )
    assert (result.initial_spg_num, result.final_spg_num) == (225, 225)
    assert result.redirected_to_symmetry is True
    assert result.n_steps == 2


def test_fixed_cell_relax_reports_max_stress() -> None:
    """Fixed-cell relaxations report the final stress as absolute magnitudes."""
    result = StandardKappaAdapter().relax(
        cu_cell(strain=0.02),
        EMT(),
        KappaSettings(ase_filter=None, max_steps=5),
        log_file=None,
    )
    assert result.max_stress is not None
    from matbench_discovery.phonons.adapters.standard import _max_stress

    (expected,) = EMT().compute([result.atoms], ["stress"])
    np.testing.assert_allclose(
        result.max_stress, _max_stress(expected["stress"]), rtol=0, atol=1e-12
    )
    assert np.abs(result.max_stress).max() > 1e-4  # strained cell is not stress-free
    # compressive stress must report its magnitude, not a signed max of ~0
    np.testing.assert_allclose(_max_stress(-3 * np.eye(3)), [3, 0], rtol=0, atol=0)


class CountingEMT:
    """EMT behind the Python model protocol, recording every call's batch size."""

    def __init__(self) -> None:
        """Start with no recorded calls."""
        self.batch_sizes: list[int] = []

    def compute(
        self, structures: list[Structure], properties: list[str]
    ) -> list[dict[str, Any]]:
        """Record the batch size, then delegate to matscikit's EMT."""
        self.batch_sizes.append(len(structures))
        return EMT().compute(structures, properties)


@pytest.mark.parametrize(
    ("batch_size", "max_atoms", "fc3_batches"),
    [(4, None, [4] * 7 + [2]), (1, 100, [12, 12, 6]), (1, 7, [1] * 30)],
)
def test_batch_limits_reach_matscikit(
    batch_size: int, max_atoms: int | None, fc3_batches: list[int]
) -> None:
    """Batches follow batch_size, or fill max_atoms_per_batch with whole supercells
    (at least one, also when a supercell alone exceeds it).
    """
    from matbench_discovery.phonons.adapters.standard import supercell_batch_size

    settings = KappaSettings(batch_size=batch_size, max_atoms_per_batch=max_atoms)
    assert supercell_batch_size(settings, 8) == fc3_batches[0]
    calc = thermal_conductivity.init_phonons(
        cu_cell(), fc2_supercell=np.eye(3), fc3_supercell=np.diag([2, 1, 1])
    )
    model = CountingEMT()
    StandardKappaAdapter().calculate_fc3(calc, model, settings)  # ty: ignore[invalid-argument-type]
    assert model.batch_sizes == fc3_batches  # 30 configurations of 8 atoms
    assert calc.phonons.fc3_forces.shape == (30, 8, 3)
