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
try:
    mp_elem_reference_entries = (
        pd.read_json(mp_elem_refs_path, typ="series").map(PDEntry.from_dict).to_dict()
    )
except FileNotFoundError:
    mp_elem_reference_entries = None


def get_e_form_per_atom(
    entry: EntryLike,
    elemental_ref_entries: dict[str, EntryLike] = None,
) -> float:
    """Get the formation energy of a composition from a list of entries and a dict
    mapping elements to reference energies.

    Args:
        entry: Entry | dict[str, float | str | Composition]: pymatgen Entry (PDEntry,
            ComputedEntry or ComputedStructureEntry) or dict with energy and composition
            keys to compute formation energy of.
        elemental_ref_entries (dict[str, Entry], optional): Must be a complete set of
            terminal (i.e. elemental) reference entries containing the lowest energy
            phase for each element present in entry. Defaults to MP elemental reference
            entries as collected on 2022-09-19 get_elemental_ref_entries(). This was
            tested to give the same formation energies as computed by MP.

    Returns:
        float: formation energy in eV/atom.
    """
    if elemental_ref_entries is None:
        if mp_elem_reference_entries is None:
            msg = f"{mp_elem_refs_path=}, pass elemental_ref_entries explicitly."
            raise FileNotFoundError(msg)
        elemental_ref_entries = mp_elem_reference_entries

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

    refs = {str(el): elemental_ref_entries[str(el)] for el in comp}

    for key, ref_entry in refs.items():
        if isinstance(ref_entry, dict):
            refs[key] = PDEntry.from_dict(ref_entry)

    form_energy = energy - sum(comp[el] * refs[str(el)].energy_per_atom for el in comp)

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
            Recall, Prevalence, Accuracy, F1, TPR, FPR, TNR, FNR, MAE, RMSE, R2.
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
        Prevalence=prevalence,
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
