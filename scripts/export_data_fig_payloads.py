"""Regenerate model-independent data-page figure payloads and element counts."""

import gzip
import json
import math
import os
from collections import Counter

import ase.io
import h5py
import pandas as pd
from ase.data import chemical_symbols
from pymatgen.core import Composition
from pymatviz.enums import Key

from matbench_discovery import MP_DIR, ROOT, SITE_FIG_DATA, figs
from matbench_discovery.data_figs import (
    build_arity_hist_payload,
    build_element_counts_payload,
    build_mp_elemental_ref_energies,
    build_mp_trj_hist_payload,
    build_route_element_counts,
    build_spacegroup_sunbursts,
    build_wbm_e_form_hist,
    build_wbm_hull_dist_hist,
)
from matbench_discovery.enums import DataFiles, MbdKey

MP_TRJ_SUMMARY_PATH = f"{MP_DIR}/2022-09-16-mp-trj-summary.json.bz2"
ROUTE_DATA_DIR = f"{ROOT}/site/src/routes/data"


def export_benchmark_element_counts() -> None:
    """Count each element once per structure, MD system, or unique diatomic pair."""
    counts = {
        "phonondb": Counter(
            symbol
            for atoms in ase.io.iread(DataFiles.phonondb_pbe_103_structures.path)
            for symbol in set(atoms.get_chemical_symbols())
        )
    }
    with h5py.File(DataFiles.dynamat_v1_0_md_trajectories.path) as reference:
        counts["dynamat"] = Counter(
            chemical_symbols[int(number)]
            for system in reference.values()
            for number in set(system["atomic_numbers"][:])
        )
    with gzip.open(DataFiles.diatomics_dft_reference.path, "rt") as file:
        # Functionals and bond-length samples share a pair; count it only once.
        pairs = {pair for curves in json.load(file).values() for pair in curves}
        counts["diatomics"] = Counter(
            symbol for pair in pairs for symbol in set(pair.split("-"))
        )
    with open(f"{ROUTE_DATA_DIR}/benchmark-element-counts.json", "w") as file:
        json.dump(counts, file, indent=2, sort_keys=True)
        file.write("\n")


def export_benchmark_summary() -> None:
    """Export public reference coverage without model predictions or withheld labels."""
    sizes = pd.read_csv(DataFiles.wbm_summary.path, usecols=["n_sites"])
    size_counts = sizes["n_sites"].value_counts().sort_index()
    with gzip.open(DataFiles.phonondb_pbe_103_kappa_no_nac.path, "rt") as file:
        phonons = json.load(file)
    kappas = [row["kappa_tot_avg"][row["temperatures"].index(300)] for row in phonons]
    if not all(math.isfinite(value) and value > 0 for value in kappas):
        raise ValueError(
            f"Expected positive finite 300 K reference conductivities: {kappas}"
        )
    systems = []
    with h5py.File(DataFiles.dynamat_v1_0_md_trajectories.path) as reference:
        for name, system in reference.items():
            numbers = system["atomic_numbers"][:]
            formula = Composition(
                Counter(chemical_symbols[int(number)] for number in numbers)
            )
            systems.append(
                {
                    "name": name,
                    "formula": formula.reduced_formula,
                    "n_atoms": len(numbers),
                    "temperature": float(system.attrs["temperature_kelvin"]),
                    "duration_ps": (system["positions"].shape[0] - 1)
                    * float(system.attrs["dt_fs"])
                    / 1000,
                }
            )
    summary = {
        "wbm_n_sites": {"x": size_counts.index.tolist(), "y": size_counts.tolist()},
        "phonondb_kappa": kappas,
        "md_systems": systems,
    }
    with open(f"{ROUTE_DATA_DIR}/benchmark-reference-summary.json", "w") as file:
        json.dump(summary, file, indent=2, allow_nan=False)
        file.write("\n")


def main() -> None:
    """Write all model-independent figure and route-local JSON payloads."""
    if not os.path.isfile(MP_TRJ_SUMMARY_PATH):
        raise FileNotFoundError(
            f"MPtrj summary cache not found at {MP_TRJ_SUMMARY_PATH}. Regenerate it "
            "from the extXYZ source with python data/mp/eda_mp_trj.py, then rerun "
            "python scripts/export_data_fig_payloads.py."
        )

    df_mp_trj = pd.read_json(MP_TRJ_SUMMARY_PATH)
    from matbench_discovery.data import df_wbm
    from matbench_discovery.energy import mp_elem_ref_entries

    df_mp = pd.read_csv(DataFiles.mp_energies.path, na_filter=False).set_index(
        Key.mat_id
    )
    wbm_spacegroups = (
        df_wbm[MbdKey.init_protostructure_spglib].str.split("_").str[2].astype(int)
    )
    mp_spacegroups = (
        df_mp[MbdKey.protostructure_spglib].str.split("_").str[2].astype(int)
    )
    payloads = {
        "hist-wbm-e-form-per-atom": build_wbm_e_form_hist(df_wbm[MbdKey.e_form_wbm]),
        "hist-wbm-hull-dist": build_wbm_hull_dist_hist(df_wbm[MbdKey.each_true]),
        "mp-elemental-ref-energies": build_mp_elemental_ref_energies(
            mp_elem_ref_entries
        ),
        "spacegroup-sunbursts": build_spacegroup_sunbursts(
            mp_spacegroups, wbm_spacegroups
        ),
        "mp-vs-mp-trj-vs-wbm-arity-hist": build_arity_hist_payload(
            df_mp[Key.formula], df_mp_trj[Key.formula], df_wbm[Key.formula]
        ),
        "mp-trj-hists": build_mp_trj_hist_payload(df_mp_trj),
    }
    for name, payload in payloads.items():
        figs.write_json_gz(f"{SITE_FIG_DATA}/{name}.json.gz", payload)

    route_counts = build_route_element_counts(df_mp, df_wbm, df_mp_trj)
    for name, counts in route_counts.items():
        counts.to_json(f"{ROUTE_DATA_DIR}/{name}.json")

    export_benchmark_element_counts()
    export_benchmark_summary()

    figs.write_json_gz(
        f"{SITE_FIG_DATA}/element-counts-mp-vs-wbm.json.gz",
        build_element_counts_payload(
            route_counts["mp-element-counts-by-occurrence"],
            df_wbm.query(MbdKey.uniq_proto)[Key.formula],
        ),
    )


if __name__ == "__main__":
    main()
