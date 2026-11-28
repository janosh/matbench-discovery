"""Conversions from WBM source entries into benchmark reference data."""

from collections.abc import Sequence
from typing import Any

import pandas as pd
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedStructureEntry
from tqdm import tqdm

from matbench_discovery.energy import mp_elemental_ref_energies
from matbench_discovery.enums import Key

# The published WBM formation-energy calculation used elemental energies to 4 decimals:
# https://github.com/janosh/matbench-discovery/blob/a7f693d1/matbench_discovery/energy.py#L74-L77
WBM_ELEMENTAL_REF_ENERGIES = {
    element: round(energy, 4) for element, energy in mp_elemental_ref_energies.items()
}


def mp2020_corrected(
    entries: Sequence[dict[str, Any]], *, verbose: bool = False
) -> pd.DataFrame:
    """Apply the current MP2020 corrections to ComputedStructureEntry dicts.

    Uses pymatgen's MaterialsProject2020Compatibility() defaults (strict_anions=
    "require_bound" since pymatgen#3803) with clean=True: corrections already stored
    on the entries are discarded.

    Args:
        entries (Sequence[dict[str, Any]]): ComputedStructureEntry dicts with
            entry_id, structure, parameters and uncorrected energy.
        verbose (bool): Whether to show a progress bar. Defaults to False.

    Returns:
        pd.DataFrame: Indexed by entry_id with columns composition (element amounts),
            uncorrected_energy (eV), correction_per_atom and energy_per_atom (eV/atom,
            corrected).

    Raises:
        ValueError: If any entry is incompatible with the MP2020 scheme.
    """
    compatibility = MaterialsProject2020Compatibility()
    entry_ids: list[str | None] = []
    rows: list[dict[str, Any]] = []
    incompatible: list[str | None] = []
    for entry_dict in tqdm(entries, desc="MP2020 corrections", disable=not verbose):
        cse = ComputedStructureEntry.from_dict(entry_dict)
        if (processed := compatibility.process_entry(cse, clean=True)) is None:
            incompatible.append(cse.entry_id)
            continue
        entry_ids.append(cse.entry_id)
        rows.append(
            {
                "composition": processed.composition.as_dict(),
                "uncorrected_energy": processed.uncorrected_energy,
                "correction_per_atom": processed.correction_per_atom,
                "energy_per_atom": processed.energy_per_atom,
            }
        )
    if incompatible:
        raise ValueError(
            f"{len(incompatible):,} entries incompatible with MP2020, first 10: "
            f"{incompatible[:10]}"
        )
    return pd.DataFrame(rows, index=pd.Index(entry_ids, name=Key.mat_id))


def unique_wbm_prototypes(summary: pd.DataFrame, mp_prototypes: pd.Series) -> pd.Series:
    """Select the lowest original WBM hull distance per novel relaxed prototype."""
    prototype = Key.protostructure_spglib
    ranked = summary.sort_values([prototype, Key.each_wbm])
    duplicate = ranked.duplicated(prototype).reindex(summary.index)
    return (~(summary[prototype].isin(mp_prototypes) | duplicate)).rename(
        Key.uniq_proto
    )
