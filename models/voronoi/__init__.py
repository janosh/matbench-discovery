import matminer.featurizers.composition as feat_comp
import matminer.featurizers.structure as feat_struct
from matminer.featurizers.base import MultipleFeaturizer

# Create the featurizer: Ward et al. use a variety of different featurizers
# https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.024104
featurizers = [
    feat_struct.SiteStatsFingerprint.from_preset("CoordinationNumber_ward-prb-2017"),
    feat_struct.StructuralHeterogeneity(),
    feat_struct.ChemicalOrdering(),
    feat_struct.MaximumPackingEfficiency(),
    feat_struct.SiteStatsFingerprint.from_preset(
        "LocalPropertyDifference_ward-prb-2017"
    ),
    feat_struct.StructureComposition(feat_comp.Stoichiometry()),
    feat_struct.StructureComposition(feat_comp.ElementProperty.from_preset("magpie")),
    feat_struct.StructureComposition(feat_comp.ValenceOrbital(props=["frac"])),
    feat_struct.StructureComposition(feat_comp.IonProperty(fast=True)),
]
featurizer = MultipleFeaturizer(featurizers)

# multiprocessing seems to be the cause of OOM errors on large structures even when
# taking only small slice of the data and launching slurm jobs with --mem 100G
featurizer.set_n_jobs(1)
