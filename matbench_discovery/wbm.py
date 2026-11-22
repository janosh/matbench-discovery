"""Conversions from WBM source entries into benchmark reference data."""

from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd
from matscikit import mp2020_corrections
from matscikit.interop import structure_from_pymatgen_dict
from pymatviz.enums import Key

from matbench_discovery.energy import mp_elemental_ref_energies
from matbench_discovery.enums import MbdKey

# The published WBM formation-energy calculation used elemental energies to 4 decimals:
# https://github.com/janosh/matbench-discovery/blob/a7f693d1/matbench_discovery/energy.py#L74-L77
WBM_ELEMENTAL_REF_ENERGIES = {
    element: round(energy, 4) for element, energy in mp_elemental_ref_energies.items()
}


def mp2020_corrected(entries: Sequence[dict[str, Any]]) -> pd.DataFrame:
    """Apply the current MP2020 corrections to ComputedStructureEntry dicts.

    Matches pymatgen's MaterialsProject2020Compatibility().process_entries(clean=True)
    (default strict_anions="require_bound") through matscikit's port: corrections
    already stored on the entries are discarded.

    Args:
        entries (Sequence[dict[str, Any]]): pymatgen-layout entry dicts with
            entry_id, structure, parameters and uncorrected energy.

    Returns:
        pd.DataFrame: Indexed by entry_id with columns composition (element amounts),
            uncorrected_energy (eV), correction_per_atom and energy_per_atom (eV/atom,
            corrected).

    Raises:
        ValueError: If any entry is incompatible with the MP2020 scheme.
    """
    structures = [structure_from_pymatgen_dict(entry["structure"]) for entry in entries]
    entry_ids = [entry["entry_id"] for entry in entries]
    results, incompatible = [], []
    for entry_id, result in zip(
        entry_ids, mp2020_corrections(entries, structures), strict=True
    ):
        if result is None:
            incompatible.append(entry_id)
        else:
            results.append(result)
    if incompatible:
        raise ValueError(
            f"{len(incompatible):,} entries incompatible with MP2020, first 10: "
            f"{incompatible[:10]}"
        )
    n_atoms = np.array([structure.composition.num_atoms for structure in structures])
    return pd.DataFrame(
        {
            "composition": [
                structure.composition.as_dict() for structure in structures
            ],
            "uncorrected_energy": [entry["energy"] for entry in entries],
            "correction_per_atom": [res.correction for res in results] / n_atoms,
            "energy_per_atom": [res.corrected_energy for res in results] / n_atoms,
        },
        index=pd.Index(entry_ids, name=Key.mat_id),
    )


def unique_wbm_prototypes(summary: pd.DataFrame, mp_prototypes: pd.Series) -> pd.Series:
    """Select the lowest original WBM hull distance per novel relaxed prototype."""
    prototype = MbdKey.protostructure_spglib
    ranked = summary.sort_values([prototype, MbdKey.each_wbm])
    duplicate = ranked.duplicated(prototype).reindex(summary.index)
    return (~(summary[prototype].isin(mp_prototypes) | duplicate)).rename(
        MbdKey.uniq_proto
    )
