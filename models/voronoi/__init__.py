import matminer.featurizers.composition as fc
import matminer.featurizers.structure as fs
from matminer.featurizers.base import MultipleFeaturizer

# Create the featurizer: Ward et al. use a variety of different featurizers
# https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.024104

composition_features = [
    # Ward+Wolverton' Magpie https://rdcu.be/c3jug
    fc.ElementProperty.from_preset("magpie"),
    # Ionic property attributes. Similar to ElementProperty.
    fc.IonProperty(fast=True),
    # Calculate norms of stoichiometric attributes.
    fc.Stoichiometry(),
    # Attributes of valence orbital shells
    fc.ValenceOrbital(props=["frac"]),
]
structure_features = [
    # How much the ordering of species in the structure differs from random
    fs.ChemicalOrdering(),
    # Maximum possible packing efficiency of this structure
    fs.MaximumPackingEfficiency(),
    # Differences in elemental properties between site and its neighboring sites
    fs.SiteStatsFingerprint.from_preset("LocalPropertyDifference_ward-prb-2017"),
    # Number of first nearest neighbors of a site.
    fs.SiteStatsFingerprint.from_preset("CoordinationNumber_ward-prb-2017"),
    # Variance in the bond lengths and atomic volumes in a structure
    fs.StructuralHeterogeneity(),
]
featurizer = MultipleFeaturizer(
    structure_features + [*map(fs.StructureComposition, composition_features)]
)


# multiprocessing seems to be the cause of OOM errors on large structures even when
# taking only small slice of the data and launching slurm jobs with --mem 100G
# Alex Dunn has been aware of this problem for a while. presumed cause: chunk of data
# (eg 50 structures) is sent to a single process, but sometimes one of those structures
# might be huge causing that process to stall. Other processes in pool can't synchronize
# at the end, effectively freezing the job
featurizer.set_n_jobs(1)
