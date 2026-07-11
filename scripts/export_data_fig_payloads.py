"""Regenerate model-independent data-page figure payloads and element counts."""

import os

import pandas as pd
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

    figs.write_json_gz(
        f"{SITE_FIG_DATA}/element-counts-mp-vs-wbm.json.gz",
        build_element_counts_payload(
            route_counts["mp-element-counts-by-occurrence"],
            df_wbm.query(MbdKey.uniq_proto)[Key.formula],
        ),
    )


if __name__ == "__main__":
    main()
