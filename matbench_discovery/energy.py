import itertools
from collections.abc import Sequence

import numpy as np
import pandas as pd
from pymatgen.analysis.phase_diagram import Entry, PDEntry
from pymatgen.core import Composition
from pymatgen.util.typing import EntryLike
from sklearn.metrics import r2_score
from tqdm import tqdm

from matbench_discovery import ROOT


def get_elemental_ref_entries(
    entries: Sequence[EntryLike], verbose: bool = True
) -> dict[str, Entry]:
    """Get the lowest energy pymatgen Entry object for each element in a list of entries.

    Args:
        entries (Sequence[Entry]): pymatgen Entries (PDEntry, ComputedEntry or
            ComputedStructureEntry) to find elemental reference entries of.
        verbose (bool, optional): Whether to show a progress bar. Defaults to False.

    Raises:
        ValueError: If some elements are missing terminal reference entries.
        ValueError: If there are more terminal entries than dimensions. Should never
            happen.

    Returns:
        dict[str, Entry]: Map from element symbol to its lowest energy entry.
    """
    entries = [PDEntry.from_dict(e) if isinstance(e, dict) else e for e in entries]
    elements = {elems for entry in entries for elems in entry.composition.elements}
    dim = len(elements)

    if verbose:
        print(f"Sorting {len(entries)} entries with {dim} dimensions...", flush=True)

    entries = sorted(entries, key=lambda e: e.composition.reduced_composition)

    elemental_ref_entries = {}
    for composition, entry_group in tqdm(
        itertools.groupby(entries, key=lambda e: e.composition.reduced_composition),
        disable=not verbose,
        desc="Finding elemental reference entries",
    ):
        min_entry = min(entry_group, key=lambda e: e.energy_per_atom)
        if composition.is_element:
            elem_symb = str(composition.elements[0])
            elemental_ref_entries[elem_symb] = min_entry

    if len(elemental_ref_entries) > dim:
        missing = elements - set(elemental_ref_entries)
        raise ValueError(f"Some terminal entries are {missing = }")
    elif len(elemental_ref_entries) < dim:
        extra = set(elemental_ref_entries) - set(elements)
        raise ValueError(
            f"There are more terminal element entries than dimensions: {extra}"
        )

    return elemental_ref_entries


# contains all MP elemental reference entries to compute formation energies
# produced by get_elemental_ref_entries() in build_phase_diagram.py
mp_elem_refs_path = f"{ROOT}/data/mp/2022-09-19-mp-elemental-reference-entries.json"
mp_elem_reference_entries = (
    pd.read_json(mp_elem_refs_path, typ="series").map(PDEntry.from_dict).to_dict()
)

# tested to agree with TRI's MP reference energies
# https://github.com/TRI-AMDD/CAMD/blob/1c965cba636531e542f4821a555b98b2d81ed034/camd/utils/data.py#L134
# fmt: off
# flake8 ignore next line
mp_elemental_ref_energies = {
    "Ne": -0.0259, "He": -0.0091, "Ar": -0.0688, "F": -1.9115, "O": -4.948, "Cl": -1.8485, "N": -8.3365, "Kr": -0.0567, "Br": -1.6369, "I": -1.524, "Xe": -0.0362, "S": -4.1364, "Se": -3.4959, "C": -9.2268, "Au": -3.2739, "W": -12.9581, "Pb": -3.7126, "Rh": -7.3643, "Pt": -6.0709, "Ru": -9.2744, "Pd": -5.1799, "Os": -11.2274, "Ir": -8.8384, "H": -3.3927, "P": -5.4133, "As": -4.6591, "Mo": -10.8456, "Te": -3.1433, "Sb": -4.129, "B": -6.6794, "Bi": -3.89, "Ge": -4.623, "Hg": -0.3037, "Sn": -4.0096, "Ag": -2.8326, "Ni": -5.7801, "Tc": -10.3606, "Si": -5.4253, "Re": -12.4445, "Cu": -4.0992, "Co": -7.1083, "Fe": -8.47, "Ga": -3.0281, "In": -2.7517, "Cd": -0.9229, "Cr": -9.653, "Zn": -1.2597, "V": -9.0839, "Tl": -2.3626, "Al": -3.7456, "Nb": -10.1013, "Be": -3.7394, "Mn": -9.162, "Ti": -7.8955, "Ta": -11.8578, "Pa": -9.5147, "U": -11.2914, "Sc": -6.3325, "Np": -12.9478, "Zr": -8.5477, "Mg": -1.6003, "Th": -7.4139, "Hf": -9.9572, "Pu": -14.2678, "Lu": -4.521, "Tm": -4.4758, "Er": -4.5677, "Ho": -4.5824, "Y": -6.4665, "Dy": -4.6068, "Gd": -14.0761, "Eu": -10.292, "Sm": -4.7186, "Nd": -4.7681, "Pr": -4.7809, "Pm": -4.7505, "Ce": -5.9331, "Yb": -1.5396, "Tb": -4.6344, "La": -4.936, "Ac": -4.1212, "Ca": -2.0056, "Li": -1.9089, "Sr": -1.6895, "Na": -1.3225, "Ba": -1.919, "Rb": -0.9805, "K": -1.1104, "Cs": -0.8954,  # noqa: E501
}
# fmt: on


def get_e_form_per_atom(
    entry: EntryLike, elemental_ref_energies: dict[str, float] = None
) -> float:
    """Get the formation energy of a composition from a list of entries and a dict
    mapping elements to reference energies.

    Args:
        entry: Entry | dict[str, float | str | Composition]: pymatgen Entry (PDEntry,
            ComputedEntry or ComputedStructureEntry) or dict with energy and composition
            keys to compute formation energy of.
        elemental_ref_entries (dict[str, float], optional): Must be a covering set (for
            entry) of terminal reference energies, i.e. eV/atom of the lowest energy
            elemental phase for each element. Defaults to MP elemental reference
            energies as collected on 2022-09-19 get_elemental_ref_entries(). This was
            tested to give the same formation energies as found in MP.

    Returns:
        float: formation energy in eV/atom.
    """
    elemental_ref_energies = elemental_ref_energies or mp_elemental_ref_energies

    if isinstance(entry, dict):
        energy = entry["energy"]
        comp = Composition(entry["composition"])  # is idempotent if already Composition
    elif isinstance(entry, Entry):
        energy = entry.energy
        comp = entry.composition
    else:
        raise TypeError(
            f"{entry=} must be Entry (or subclass like ComputedEntry) or dict"
        )

    e_refs = {str(el): elemental_ref_energies[str(el)] for el in comp}

    for key, ref_entry in e_refs.items():
        if isinstance(ref_entry, dict):
            e_refs[key] = PDEntry.from_dict(ref_entry)

    form_energy = energy - sum(comp[el] * e_refs[str(el)] for el in comp)

    return form_energy / comp.num_atoms


def classify_stable(
    e_above_hull_true: pd.Series,
    e_above_hull_pred: pd.Series,
    stability_threshold: float = 0,
) -> tuple[pd.Series, pd.Series, pd.Series, pd.Series]:
    """Classify model stability predictions as true/false positive/negatives (usually
    w.r.t DFT-ground truth labels). All energies are assumed to be in eV/atom
    (but shouldn't really matter as long as they're consistent).

    Args:
        e_above_hull_true (pd.Series): Ground truth energy above convex hull values.
        e_above_hull_pred (pd.Series): Model predicted energy above convex hull values.
        stability_threshold (float, optional): Maximum energy above convex hull for a
            material to still be considered stable. Usually 0, 0.05 or 0.1. Defaults to
            0, meaning a material has to be directly on the hull to be called stable.
            Negative values mean a material has to pull the known hull down by that
            amount to count as stable. Few materials lie below the known hull, so only
            negative values very close to 0 make sense.

    Returns:
        tuple[TP, FN, FP, TN]: Indices as pd.Series for true positives,
            false negatives, false positives and true negatives (in this order).
    """
    actual_pos = e_above_hull_true <= stability_threshold
    actual_neg = e_above_hull_true > stability_threshold
    model_pos = e_above_hull_pred <= stability_threshold
    model_neg = e_above_hull_pred > stability_threshold

    true_pos = actual_pos & model_pos
    false_neg = actual_pos & model_neg
    false_pos = actual_neg & model_pos
    true_neg = actual_neg & model_neg

    return true_pos, false_neg, false_pos, true_neg


def stable_metrics(
    true: Sequence[float], pred: Sequence[float], stability_threshold: float = 0
) -> dict[str, float]:
    """
    Get a dictionary of stability prediction metrics. Mostly binary classification
    metrics, but also MAE, RMSE and R2.

    Args:
        true (list[float]): true energy values
        pred (list[float]): predicted energy values
        stability_threshold (float): Where to place stability threshold relative to
            convex hull in eV/atom, usually 0 or 0.1 eV. Defaults to 0.

    Note: Could be replaced by sklearn.metrics.classification_report() which takes
        binary labels. I.e. classification_report(true > 0, pred > 0, output_dict=True)
        should give equivalent results.

    Returns:
        dict[str, float]: dictionary of classification metrics with keys DAF, Precision,
            Recall, Accuracy, F1, TPR, FPR, TNR, FNR, MAE, RMSE, R2.
    """
    true_pos, false_neg, false_pos, true_neg = classify_stable(
        true, pred, stability_threshold
    )

    n_true_pos, n_false_pos, n_true_neg, n_false_neg = map(
        sum, (true_pos, false_pos, true_neg, false_neg)
    )

    n_total_pos = n_true_pos + n_false_neg
    prevalence = n_total_pos / len(true)  # null rate
    precision = n_true_pos / (n_true_pos + n_false_pos)
    recall = n_true_pos / n_total_pos

    is_nan = np.isnan(true) | np.isnan(pred)
    true, pred = np.array(true)[~is_nan], np.array(pred)[~is_nan]

    return dict(
        DAF=precision / prevalence,
        Precision=precision,
        Recall=recall,
        Accuracy=(n_true_pos + n_true_neg) / len(true),
        F1=2 * (precision * recall) / (precision + recall),
        TPR=n_true_pos / (n_true_pos + n_false_neg),
        FPR=n_false_pos / (n_true_neg + n_false_pos),
        TNR=n_true_neg / (n_true_neg + n_false_pos),
        FNR=n_false_neg / (n_true_pos + n_false_neg),
        MAE=np.abs(true - pred).mean(),
        RMSE=((true - pred) ** 2).mean() ** 0.5,
        R2=r2_score(true, pred),
    )
