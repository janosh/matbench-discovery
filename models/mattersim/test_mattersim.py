"""
Get MatterSim predictions on WBM dataset

Please note this script does not run because we removed
the codes that may disclose any usage information of MatterSim.
"""

from __future__ import annotations

import bz2
import json
from typing import TYPE_CHECKING

import pandas as pd
from ase.constraints import ExpCellFilter
from ase.optimize import FIRE
from pymatgen.core import Structure
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.data import DataFiles
from matbench_discovery.energy import get_e_form_per_atom

if TYPE_CHECKING:
    from typing import Literal

    import ase
    from ase.calculators.singlepoint import SinglePointCalculator

__author__ = "Han Yang"
__date__ = "2024-06-19"


def dummy_mattersim_calculator(
    backbone: Literal["m3gnet", "graphormer"] = "m3gnet",
) -> SinglePointCalculator:
    """
    This is a dummy function that makes a MatterSim calculator
    """
    print(backbone)
    calc: SinglePointCalculator = None
    return calc


def convert_wbm_to_atoms_list() -> list[ase.Atoms]:
    """
    This function reads the initial structures from the
    WBM dataset and converts them to ASE atoms objects
    """
    with bz2.open(
        "/path/to/initial/wbm/structures/2022-10-19-wbm-init-structs.json.bz2"
    ) as fh:
        data_initial = json.loads(fh.read().decode("utf-8"))

    atoms_list_initial = []
    for key, structure_dict in tqdm(data_initial["initial_structure"].items()):
        structure = Structure.from_dict(structure_dict)
        atoms = AseAtomsAdaptor.get_atoms(structure)
        atoms.info["material_id"] = data_initial["material_id"][key]
        atoms_list_initial.append(atoms)

    return atoms_list_initial


def relax_atoms_list(
    atoms_list: list[ase.Atoms],
    fmax: float = 0.01,
    steps: int = 500,
    backbone: Literal["m3gnet", "graphormer"] = "m3gnet",
) -> list[ase.Atoms]:
    """
    This function relax the atoms.
    """
    calc_m3gnet = dummy_mattersim_calculator(backbone="m3gnet")

    relaxed_atoms_list = []

    for atoms in atoms_list:
        atoms.set_calculator(calc_m3gnet)
        ecf = ExpCellFilter(atoms)
        opt = FIRE(ecf)
        opt.run(fmax=fmax, steps=steps)
        if opt.get_number_of_steps() == steps:
            atoms.info["converged"] = False
        else:
            atoms.info["converged"] = True

        if backbone == "graphormer":
            # Please note that we only re-calculate the
            # energy in the case of MatterSim(graphormer).
            # The structure relaxation is always done with MatterSim(m3gnet).
            calc = dummy_mattersim_calculator(backbone="graphormer")
            atoms.set_calculator(calc)

        relaxed_atoms_list.append(atoms)

    return relaxed_atoms_list


# %%
def parse_relaxed_atoms_list_as_df(
    atoms_list: list[ase.Atoms], *, keep_unconverged: bool = True
) -> pd.DataFrame:
    e_form_col = "e_form_per_atom_mattersim"

    wbm_cse_paths = DataFiles.wbm_computed_structure_entries.path
    df_cse = pd.read_json(wbm_cse_paths).set_index(Key.mat_id)

    df_cse[Key.cse] = [
        ComputedStructureEntry.from_dict(dct) for dct in tqdm(df_cse[Key.cse])
    ]

    print(f"Found {len(df_cse):,} CSEs in {wbm_cse_paths=}")
    print(f"Found {len(atoms_list):,} relaxed structures")

    def parse_single_atoms(atoms: ase.Atoms) -> tuple[str, bool, float, float, float]:
        structure = AseAtomsAdaptor.get_structure(atoms)
        energy = atoms.get_potential_energy()
        mat_id = atoms.info["material_id"]
        converged = atoms.info["converged"]

        cse = df_cse.loc[mat_id, Key.cse]
        cse._energy = energy  # noqa: SLF001
        cse._structure = structure  # noqa: SLF001

        processed = MaterialsProject2020Compatibility(check_potcar=False).process_entry(
            cse, verbose=False, clean=True
        )
        corrected_energy = processed.energy if processed is not None else energy
        formation_energy = (
            get_e_form_per_atom(processed)
            if processed is not None
            else get_e_form_per_atom(cse)
        )

        return mat_id, converged, formation_energy, energy, corrected_energy

    mat_id_list, converged_list, e_form_list = [], [], []
    energy_list, corrected_energy_list = [], []

    for atoms in tqdm(atoms_list, "Processing relaxed structures"):
        mat_id, converged, formation_energy, energy, corrected_energy = (
            parse_single_atoms(atoms)
        )
        if not keep_unconverged and not converged:
            continue
        mat_id_list += [mat_id]
        converged_list += [converged]
        e_form_list += [formation_energy]
        energy_list += [energy]
        corrected_energy_list += [corrected_energy]

    return pd.DataFrame(
        {
            Key.mat_id: mat_id_list,
            "converged": converged_list,
            e_form_col: e_form_list,
            "mattersim_energy": energy_list,
            "corrected_energy": corrected_energy_list,
        }
    )


if __name__ == "__main__":
    init_wbm_atoms_list = convert_wbm_to_atoms_list()
    relaxed_wbm_atoms_list = relax_atoms_list(
        init_wbm_atoms_list, backbone="graphormer"
    )
    parse_relaxed_atoms_list_as_df(relaxed_wbm_atoms_list).to_csv(
        f"{today}-mattersim-wbm-IS2RE.csv.gz"
    )
