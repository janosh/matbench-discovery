"""Recreate the featurizer used by Ward et al. in
https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.024104.
"""

import matminer.featurizers.composition as fc
import matminer.featurizers.structure as fs
from matminer.featurizers.base import MultipleFeaturizer

composition_features = [
    # Ward+Wolverton' Magpie https://rdcu.be/c3jug
    fc.ElementProperty.from_preset("magpie"),
    fc.IonProperty(fast=True),
    fc.Stoichiometry(),
    # Attributes of valence orbital shells
    fc.ValenceOrbital(props=["frac"]),
]
structure_features = [
    # How much the ordering of species in the structure differs from random
    fs.ChemicalOrdering(),
    fs.MaximumPackingEfficiency(),
    # Differences in elemental properties between site and its neighboring sites
    fs.SiteStatsFingerprint.from_preset("LocalPropertyDifference_ward-prb-2017"),
    # Number of first nearest neighbors of a site.
    fs.SiteStatsFingerprint.from_preset("CoordinationNumber_ward-prb-2017"),
    # Variance in the bond lengths and atomic volumes in a structure
    fs.StructuralHeterogeneity(),
]
featurizer = MultipleFeaturizer(
    [*structure_features, *map(fs.StructureComposition, composition_features)]
)


# multiprocessing OOMs on large structures even on small data slices with --mem 100G
# (long known to Alex Dunn). Presumed cause: a chunk (eg 50 structures) goes to one
# process, but a single huge structure stalls it, so the pool never synchronizes at the
# end and the job freezes.
featurizer.set_n_jobs(1)
