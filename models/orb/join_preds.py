import glob

import pandas as pd
import typer
from pymatgen.core import Structure
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery.energy import get_e_form_per_atom
from matbench_discovery.enums import MbdKey

app = typer.Typer(pretty_exceptions_enable=False, no_args_is_help=True)
FORMATION_ENERGY_COL = "e_form_per_atom_orb"
STRUCT_COL = "orb_structure"


@app.command()
def main(
    predictions_dir: str,
    glob_pattern: str = "*shard*.json.gz",
    correct_energies: bool = True,  # noqa: FBT001,FBT002
) -> None:
    """Join ORB predictions with WBM data.

    In addition, this script applies the MaterialsProject2020Compatibility
    energy corrections to the ORB predictions and computes the corrected
    formation energies.

    This script produces 4 files:
    - {predictions_dir}/{prefix}.csv.gz (all predictions)
    - {predictions_dir}/{prefix}-no-bad.csv.gz
      (all predictions except ones with large errors)
    - {predictions_dir}/{prefix}.json.gz (all predictions as JSON)
    - {predictions_dir}/{prefix}-bad.csv (predictions with large errors)
    """
    file_paths = sorted(glob.glob(f"{predictions_dir}/{glob_pattern}"))

    print(f"Found {len(file_paths):,} files for {glob_pattern = }")
    dfs: dict[str, pd.DataFrame] = {}

    for file_path in tqdm(file_paths):
        if file_path in dfs:
            continue
        df_i = pd.read_json(file_path).set_index(Key.mat_id)
        # drop trajectory to save memory, if present
        dfs[file_path] = df_i.drop(columns="orb_trajectory", errors="ignore")

    df_orb = pd.concat(dfs.values()).round(4)

    # This is inside the script because accessing the variables causes a download
    # to be triggered if they are not present, meaning it's better to only load them
    # if the script is actually going to be run.
    from matbench_discovery.data import DataFiles, as_dict_handler, df_wbm

    if correct_energies:
        df_cse = pd.read_json(DataFiles.wbm_computed_structure_entries.path).set_index(
            Key.mat_id
        )

        df_cse[Key.cse] = [
            ComputedStructureEntry.from_dict(dct)
            for dct in tqdm(df_cse[Key.cse], desc="Loading CSEs")
        ]

        # transfer predicted energies and relaxed structures WBM CSEs since
        # MP2020 energy corrections applied below are structure-dependent
        # (for oxides and sulfides)
        cse: ComputedStructureEntry
        for row in tqdm(
            df_orb.itertuples(), total=len(df_orb), desc="ML energies to CSEs"
        ):
            mat_id, struct_dict, orb_energy, *_ = row
            mlip_struct = Structure.from_dict(struct_dict)
            df_orb.loc[mat_id, STRUCT_COL] = mlip_struct
            cse = df_cse.loc[mat_id, Key.cse]
            # cse._energy is the uncorrected energy
            cse._energy = orb_energy  # noqa: SLF001
            cse._structure = mlip_struct  # noqa: SLF001
            df_orb.loc[mat_id, Key.cse] = cse

        # apply energy corrections inplace
        processed = MaterialsProject2020Compatibility().process_entries(
            df_orb[Key.cse], verbose=True, clean=True
        )
        if len(processed) != len(df_orb):
            raise ValueError(
                f"not all entries processed: {len(processed)=} {len(df_orb)=}"
            )

        # compute corrected formation energies
        df_orb[Key.formula] = df_wbm[Key.formula]
        df_orb[FORMATION_ENERGY_COL] = [
            get_e_form_per_atom(dict(energy=cse.energy, composition=formula))
            for formula, cse in tqdm(
                df_orb.set_index(Key.formula)[Key.cse].items(), total=len(df_orb)
            )
        ]

    else:
        df_orb[Key.formula] = df_wbm[Key.formula]
        df_orb[FORMATION_ENERGY_COL] = [
            get_e_form_per_atom(dict(energy=energy, composition=formula))
            for formula, energy in tqdm(
                df_orb.set_index(Key.formula)["orb_energy"].items(), total=len(df_orb)
            )
        ]

    df_wbm[[*df_orb]] = df_orb

    bad_mask = (df_wbm[FORMATION_ENERGY_COL] - df_wbm[MbdKey.e_form_dft]) < -5
    print(f"{sum(bad_mask)=}")

    # e.g orbFF-v1-2024-07-11-shard-001.json.gz -> orbFF-v1-2024-07-11
    first_filename = file_paths[0].rsplit("/", 1)[1]
    prefix = first_filename.split("-shard")[0]
    if not correct_energies:
        prefix += "-no-corr"
    out_path = f"{predictions_dir}/{prefix}"

    df_orb = df_orb.round(4)
    df_orb.select_dtypes("number").to_csv(f"{out_path}.csv.gz")
    df_orb[~bad_mask].select_dtypes("number").to_csv(f"{out_path}-no-bad.csv.gz")
    df_orb.reset_index().to_json(f"{out_path}.json.gz", default_handler=as_dict_handler)

    df_bad = df_orb[bad_mask].drop(columns=[Key.cse, STRUCT_COL], errors="ignore")
    df_bad[MbdKey.e_form_dft] = df_wbm[MbdKey.e_form_dft]
    df_bad.to_csv(f"{out_path}-bad.csv")


if __name__ == "__main__":
    app()
