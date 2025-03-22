"""Functions to analyze symmetry of sets of structures."""

import pandas as pd
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery.enums import MbdKey


def get_sym_info_from_structs(
    structures: dict[str, Structure],
    *,
    pbar: bool | dict[str, str | float | bool] = True,
    symprec: float = 1e-2,
    angle_tolerance: float | None = None,
) -> pd.DataFrame:
    """Compile DataFrame of high-level symmetry information for a dictionary of
    structures.

    Args:
        structures (dict[str, Structure | Atoms]): Map of material IDs to pymatgen
            Structures or ASE Atoms objects
        pbar (bool | dict[str, str | float | bool], optional): Whether to show progress
            bar. Defaults to True.
        symprec (float, optional): Symmetry precision of moyopy. Defaults to 1e-2.
        angle_tolerance (float, optional): Angle tolerance of moyopy (in radians unlike
            spglib which uses degrees!). Defaults to None.

    Returns:
        pd.DataFrame: DataFrame containing symmetry information for each structure
    """
    import moyopy
    from moyopy.interface import MoyoAdapter

    results: dict[str, dict[str, str | int | list[str]]] = {}
    iterator = structures.items()
    if pbar:
        pbar_kwargs = pbar if isinstance(pbar, dict) else {}
        pbar_kwargs.setdefault("desc", "Analyzing symmetry")
        iterator = tqdm(iterator, total=len(structures), **pbar_kwargs)

    for struct_key, struct in iterator:
        moyo_cell = MoyoAdapter.from_py_obj(struct)

        sym_data = moyopy.MoyoDataset(
            moyo_cell, symprec=symprec, angle_tolerance=angle_tolerance
        )

        sym_ops = sym_data.operations
        hall_symbol_entry = moyopy.HallSymbolEntry(hall_number=sym_data.hall_number)

        sym_info = {
            Key.spg_num: sym_data.number,
            Key.hall_num: sym_data.hall_number,
            MbdKey.international_spg_name: sym_data.site_symmetry_symbols,
            Key.wyckoff_symbols: sym_data.wyckoffs,
            Key.n_sym_ops: sym_ops.num_operations,
            Key.n_rot_syms: len(sym_ops.rotations),
            Key.n_trans_syms: len(sym_ops.translations),
            Key.hall_symbol: hall_symbol_entry.hm_short,
            Key.hall_num: hall_symbol_entry.hall_number,
        }
        results[struct_key] = sym_info | dict(
            symprec=symprec, angle_tolerance=angle_tolerance
        )

    df_sym = pd.DataFrame(results).T
    df_sym.index.name = Key.mat_id
    return df_sym


def calc_structure_distances(
    df_result: pd.DataFrame,
    pred_structs: dict[str, Structure],
    ref_structs: dict[str, Structure],
    *,
    pbar: bool | dict[str, str | float | bool] = True,
) -> pd.DataFrame:
    """Calculate distance metrics between prediction and reference structures.

    Args:
        df_result (pd.DataFrame): DataFrame to which distance metrics will be added.
        pred_structs (dict[str, Structure]): Map material IDs to prediction structures.
        ref_structs (dict[str, Structure]): Map material IDs to reference structures.
        pbar (bool | dict[str, str | float | bool], optional): Progress bar config.
            Defaults to True.

    Returns:
        pd.DataFrame: DataFrame with added columns for structure distance metrics
    """
    structure_matcher = StructureMatcher(stol=1.0, scale=False)
    ref_ids, pred_ids = set(ref_structs), set(pred_structs)
    shared_ids = ref_ids & pred_ids

    # Initialize RMSD column
    df_result[MbdKey.structure_rmsd_vs_dft] = None
    df_result[Key.max_pair_dist] = None

    if len(shared_ids) == 0:
        print(
            f"⚠️ Warning: No shared IDs between {len(pred_ids)} model and "
            f"{len(ref_ids)} reference structures."
        )
        print(f"First 5 model IDs: {list(pred_ids)[:5]}")
        print(f"First 5 reference IDs: {list(ref_ids)[:5]}")
        return df_result  # Return DataFrame with NaN values for RMSD

    if pbar:
        pbar_kwargs = dict(leave=False, desc="Calculating RMSD") | (
            {} if pbar is True else pbar
        )
        shared_ids = tqdm(shared_ids, **pbar_kwargs)

    for mat_id in shared_ids:
        if mat_id in df_result.index:
            rmsd, max_dist = structure_matcher.get_rms_dist(
                pred_structs[mat_id], ref_structs[mat_id]
            ) or (None, None)
            df_result.loc[mat_id, MbdKey.structure_rmsd_vs_dft] = rmsd
            df_result.loc[mat_id, Key.max_pair_dist] = max_dist

    return df_result
