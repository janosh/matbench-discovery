"""Formation and elemental reference energies from composition-energy entries."""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from matscikit import Composition, PhaseDiagram

from matbench_discovery.enums import DataFiles


@dataclass(frozen=True, init=False)
class Entry:
    """A composition with its total (corrected) energy in eV."""

    composition: Composition
    energy: float
    entry_id: str | None

    def __init__(
        self,
        composition: Composition | str | Mapping[str, float],
        energy: float,
        entry_id: str | None = None,
    ) -> None:
        """Coerce a formula string or element->amount mapping to a Composition."""
        if not isinstance(composition, Composition):
            composition = Composition(composition)
        object.__setattr__(self, "composition", composition)
        object.__setattr__(self, "energy", float(energy))
        object.__setattr__(self, "entry_id", entry_id)

    @property
    def energy_per_atom(self) -> float:
        """Total energy divided by the number of atoms."""
        return self.energy / self.composition.num_atoms

    @classmethod
    def from_dict(cls, dct: Mapping[str, Any]) -> Entry:
        """Build from a serialized entry (pymatgen PDEntry/ComputedEntry layout).

        ComputedEntry dicts store the uncorrected energy plus the total ``correction``
        (the sum of their energy_adjustments), so the corrected energy is
        ``energy + correction`` exactly as pymatgen's ``ComputedEntry.energy``.
        """
        energy = float(dct["energy"]) + float(dct.get("correction", 0.0))
        return cls(dct["composition"], energy, dct.get("entry_id"))


# contains all MP elemental reference entries to compute formation energies
# produced by scripts/data/build_mp_reference_entries.py
mp_elem_ref_entries = (
    pd.read_json(DataFiles.mp_elemental_ref_entries.path, typ="series")
    .map(Entry.from_dict)
    .to_dict()
)

# tested to agree with TRI's MP reference energies
# https://github.com/TRI-AMDD/CAMD/blob/1c965cba636/camd/utils/data.py#L134
mp_elemental_ref_energies = {
    elem: entry.energy_per_atom for elem, entry in mp_elem_ref_entries.items()
}


def formation_energies(
    compositions: Sequence[str | Composition | Mapping[str, float]],
    total_energies: Sequence[float] | np.ndarray,
    ref_energies: Mapping[str, float],
) -> np.ndarray:
    """Energies per atom (eV/atom) relative to elemental reference energies per atom.

    With elemental ground-state references (e.g. mp_elemental_ref_energies) these
    are formation energies, with isolated-atom references cohesive energies.

    Args:
        compositions (Sequence[str | Composition | Mapping[str, float]]): Formula
            strings, Compositions or element->amount mappings.
        total_energies (Sequence[float] | np.ndarray): Total energies (eV), one per
            composition.
        ref_energies (Mapping[str, float]): Reference energy per atom of each element.

    Raises:
        ValueError: If a composition has an element without a reference energy.

    Returns:
        np.ndarray: Energy per atom relative to the references, one per composition.
    """
    n_atoms = [Composition(composition).num_atoms for composition in compositions]
    elemental_refs = PhaseDiagram(list(ref_energies), list(ref_energies.values()))
    return np.asarray(
        elemental_refs.formation_energy(
            list(compositions), np.divide(total_energies, n_atoms).tolist()
        )
    )
