"""Check electronic-state validation and ASE caching without optional GPU packages."""

import importlib.util
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator


@pytest.fixture
def calculator(monkeypatch: pytest.MonkeyPatch) -> Calculator:
    """Load the real wrapper with a stand-in for the optional foundation runtime."""
    upstream = ModuleType("mace_soyo")
    monkeypatch.setattr(upstream, "MACESoyoCalculator", Calculator, raising=False)
    monkeypatch.setitem(sys.modules, "mace_soyo", upstream)
    monkeypatch.setitem(sys.modules, "torch", ModuleType("torch"))
    path = (
        Path(__file__).resolve().parents[2]
        / "matbench_discovery/mace_soyo_hyd/calculator.py"
    )
    spec = importlib.util.spec_from_file_location("hyd_input_validation", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    model = SimpleNamespace(_conditioned=False, coefficients=np.zeros((4, 4)))
    monkeypatch.setattr(module, "MaceSoyoHydCalculator", lambda *_args: model)
    calc = module.CanonicalBatchedDimerCalculator("unused", "unused", "cpu")
    monkeypatch.setattr(
        calc,
        "predict",
        lambda group: [
            {"energy": -1.25, "free_energy": -1.25, "forces": np.zeros((2, 3))}
            for _ in group
        ],
    )
    return calc


@pytest.mark.parametrize("state", [{"charge": 0, "spin": 0}, {"charge": 1, "spin": 2}])
def test_unsupported_state_rejected_before_geometry_cache(
    calculator: Calculator, state: dict[str, int]
) -> None:
    """Cached neutral geometry must never supply results for unsupported metadata."""
    atoms = Atoms("HF", positions=[[0, 0, 0], [1, 0, 0]], calculator=calculator)
    assert atoms.get_potential_energy() == -1.25
    assert calculator.cached
    atoms.info.update(state)
    with pytest.raises(ValueError, match="package has no charge/spin inputs"):
        atoms.get_potential_energy()
    atoms.info.clear()
    assert atoms.get_potential_energy() == -1.25


@pytest.mark.parametrize("name", ["charge", "spin"])
@pytest.mark.parametrize("previous, current", [(None, 0), (0, 1), (0, None)])
def test_electronic_state_changes_invalidate_ase_results(
    calculator: Calculator, name: str, previous: int | None, current: int | None
) -> None:
    """ASE geometry equality alone is insufficient for electronic-state changes."""
    atoms = Atoms("HF", positions=[[0, 0, 0], [1, 0, 0]], info={name: previous})
    calculator.atoms = atoms.copy()
    assert calculator.check_state(atoms) == []
    atoms.info[name] = current
    assert name in calculator.check_state(atoms)


@pytest.mark.parametrize("state", [{}, {"charge": 0}, {"spin": 1}])
def test_neutral_and_partial_metadata_keep_upstream_contract(
    calculator: Calculator, state: dict[str, int]
) -> None:
    """Unspecified conditioning and partial metadata retain upstream behavior."""
    atoms = Atoms("HF", positions=[[0, 0, 0], [1, 0, 0]], info=state)
    calculator.calculate(atoms)
    assert calculator.results["energy"] == -1.25
    np.testing.assert_array_equal(calculator.results["forces"], np.zeros((2, 3)))
    assert calculator.check_state(atoms) == []


def test_direct_batch_rejects_before_input_conversion(
    calculator: Calculator, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Public batch calls validate every member before any model conversion or call."""
    neutral = Atoms("HF", positions=[[0, 0, 0], [1, 0, 0]])
    unsupported = neutral.copy()
    unsupported.info.update(charge=0, spin=0)
    converted = []
    monkeypatch.setattr(
        calculator.model, "_atoms_to_inputs", converted.append, raising=False
    )
    with pytest.raises(ValueError, match="package has no charge/spin inputs"):
        calculator.predict_base([neutral, unsupported])
    assert converted == []
