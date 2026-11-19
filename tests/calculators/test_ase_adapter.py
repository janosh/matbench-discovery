"""Pin matscikit.interop.AseCalculatorModel as a transparent wrapper of ASE calculators.

Tolerances: structure_from_ase stores fractional coordinates and structure_to_ase
rebuilds Cartesian ones, so the adapter's atoms differ from the input by a few ulp
(measured <= 1.8e-15 Å on |x| <= 10.6 Å, i.e. <= 1.6 eps max|x|; bound 4 eps). Fed the
round-tripped atoms, the adapter must match the bare calculator bit for bit. Against
the original atoms, results shift by dE ~ F dx and dF ~ H dx. Measured maxima over
these structures (EMT / MACE-MP-0 in float64): 2.0e-15 / 1.2e-15 eV/atom,
7.3e-14 / 7.8e-14 eV/Å and 4.5e-16 / 4.1e-18 eV/Å³. The bounds below sit 13-50x above
them, still >= 10x tighter than any physically meaningful difference.
"""

from collections.abc import Callable
from typing import Any

import numpy as np
import pytest
from matscikit.interop import AseCalculatorModel, structure_from_ase, structure_to_ase

ase = pytest.importorskip("ase")
from ase import Atoms  # noqa: E402
from ase.build import bulk, fcc111, molecule  # noqa: E402
from ase.calculators.calculator import Calculator  # noqa: E402
from ase.calculators.emt import EMT  # noqa: E402
from ase.geometry import cellpar_to_cell  # noqa: E402

F64_EPS = np.finfo(np.float64).eps
POSITION_ULPS = 4  # round-trip bound in eps * max(|x|, 1)
ENERGY_ATOL = 1e-13  # eV/atom
FORCE_ATOL = 1e-12  # eV/Å
STRESS_ATOL = 1e-14  # eV/Å³
EMT_ELEMENTS = ("Al", "Ni", "Cu", "Pd", "Ag", "Pt", "Au")
PROPERTIES = ("energy", "forces", "stress")


def _random_triclinic(rng: np.random.Generator, n_atoms: int) -> Atoms:
    """Random triclinic cell (angles 65-115°, edges >= 2.5 Å, ~18 Å³/atom) of
    EMT elements placed sequentially at least 2 Å apart under minimum image.
    """
    lengths = np.maximum(rng.uniform(0.85, 1.2, 3) * (18 * n_atoms) ** (1 / 3), 2.5)
    cell = cellpar_to_cell([*lengths, *rng.uniform(65, 115, 3)])
    atoms = Atoms(cell=cell, pbc=True)
    while len(atoms) < n_atoms:
        trial = atoms + Atoms(
            [rng.choice(EMT_ELEMENTS)], scaled_positions=rng.random((1, 3)), cell=cell
        )
        others = range(len(trial) - 1)
        if not others or trial.get_distances(-1, others, mic=True).min() > 2.0:
            atoms = trial
    return atoms


def _make_structures() -> dict[str, Atoms]:
    """Seeded triclinic cells plus a 1-atom cell with many periodic images inside
    the EMT cutoff, a slab with mixed pbc and an isolated molecule in a box.
    """
    rng = np.random.default_rng(seed=0)
    structs = {
        f"triclinic-{idx}": _random_triclinic(rng, int(n_atoms))
        for idx, n_atoms in enumerate(rng.integers(1, 9, size=30))
    }
    structs["cu-primitive-many-images"] = bulk("Cu", "fcc", a=3.6)
    slab = fcc111("Pt", size=(2, 2, 3), vacuum=6.0)
    slab.pbc = (True, True, False)
    slab.rattle(stdev=0.05, seed=0)
    structs["pt111-slab-mixed-pbc"] = slab
    ethanol = molecule("CH3CH2OH")
    ethanol.center(vacuum=5.0)  # pbc stays False, box only sets the stress volume
    structs["ethanol-isolated"] = ethanol
    return structs


STRUCTURES = _make_structures()
# MACE on CPU costs ~60 ms per structure: every non-random case plus 10 random cells
MACE_KEYS = list(dict.fromkeys([*list(STRUCTURES)[:30:3], *list(STRUCTURES)[30:]]))
CASES = [(key, "emt") for key in STRUCTURES] + [(key, "mace") for key in MACE_KEYS]


def _direct(atoms: Atoms, calc: Calculator) -> dict[str, Any]:
    """Energy, forces and 3x3 stress from the calculator attached to ``atoms``."""
    atoms = atoms.copy()
    atoms.calc = calc
    return {
        "energy": atoms.get_potential_energy(),
        "forces": atoms.get_forces(),
        "stress": atoms.get_stress(voigt=False),
    }


@pytest.fixture(scope="module")
def mace_calc() -> Calculator:
    """MACE-MP-0 medium on CPU in float64 (skips without mace-torch)."""
    mace_calcs = pytest.importorskip("mace.calculators")
    return mace_calcs.mace_mp(model="medium", device="cpu", default_dtype="float64")


@pytest.mark.parametrize(("struct_key", "calc_name"), CASES)
def test_ase_calculator_model_matches_direct_calculator(
    struct_key: str, calc_name: str, request: pytest.FixtureRequest
) -> None:
    """Adapter equals the bare ASE calculator up to the position round trip."""
    atoms = STRUCTURES[struct_key]
    calc = EMT() if calc_name == "emt" else request.getfixturevalue("mace_calc")

    structure = structure_from_ase(atoms)
    round_tripped = structure_to_ase(structure)
    pos_tol = POSITION_ULPS * F64_EPS * max(1.0, np.abs(atoms.positions).max())
    assert np.abs(round_tripped.positions - atoms.positions).max() <= pos_tol
    np.testing.assert_array_equal(round_tripped.cell[:], atoms.cell[:])
    np.testing.assert_array_equal(round_tripped.pbc, atoms.pbc)
    np.testing.assert_array_equal(round_tripped.numbers, atoms.numbers)

    (result,) = AseCalculatorModel(calc).compute([structure], PROPERTIES)
    assert set(result) == set(PROPERTIES)
    assert result["stress"].shape == (3, 3)

    # same input atoms: the adapter adds nothing, bit for bit
    same_input = _direct(round_tripped, calc)
    assert result["energy"] == same_input["energy"]
    np.testing.assert_array_equal(result["forces"], same_input["forces"])
    np.testing.assert_array_equal(result["stress"], same_input["stress"])

    # original atoms: differences only from the few-ulp position round trip
    direct = _direct(atoms, calc)
    np.testing.assert_allclose(
        result["energy"] / len(atoms),
        direct["energy"] / len(atoms),
        rtol=0,
        atol=ENERGY_ATOL,
    )
    np.testing.assert_allclose(
        result["forces"], direct["forces"], rtol=0, atol=FORCE_ATOL
    )
    np.testing.assert_allclose(
        result["stress"], direct["stress"], rtol=0, atol=STRESS_ATOL
    )


def test_ase_calculator_model_keeps_info_drops_arrays() -> None:
    """atoms.info reaches the calculator, per-atom arrays don't."""
    # harmless while no matbench calculator reads magmoms, charges or tags
    seen: list[Atoms] = []

    class SpyEMT(EMT):
        """EMT recording a copy of the atoms each calculation sees."""

        def calculate(
            self, atoms: Atoms | None = None, *args: object, **kwargs: object
        ) -> None:
            """Record ``atoms`` then defer to EMT."""
            assert atoms is not None
            seen.append(atoms.copy())
            super().calculate(atoms, *args, **kwargs)

    atoms = bulk("Cu", "fcc", a=3.6, cubic=True)
    info = {"material_id": "test-1", "n_steps": 3, "temperature_k": 300.5}
    atoms.info.update(info)
    n_atoms = len(atoms)
    dropped_arrays: dict[str, Callable[[np.ndarray], None]] = {
        "initial_magmoms": atoms.set_initial_magnetic_moments,
        "initial_charges": atoms.set_initial_charges,
        "tags": atoms.set_tags,
    }
    for name, setter in dropped_arrays.items():
        setter(np.arange(1, n_atoms + 1, dtype=int if name == "tags" else float))
        assert name in atoms.arrays

    AseCalculatorModel(SpyEMT()).compute([structure_from_ase(atoms)], ["energy"])

    (calc_atoms,) = seen
    assert {key: calc_atoms.info[key] for key in info} == info
    assert set(calc_atoms.arrays) == {"numbers", "positions"}
    assert not set(dropped_arrays) & set(calc_atoms.arrays)
