# Kairos Materials

Prophet is a family of atomistic foundation models for materials simulation and discovery, trained on large-scale first-principles data to model energies, forces and stresses across diverse chemical systems, with extensions that explicitly incorporate spin as a fundamental physical degree of freedom. Prophet-OAME-MBD is a 62.3M-parameter model pretrained on Omat24 and ELEMENTA, including the ELEMENTA Vib expansion, and aligned on MPtrj and sAlex to the Materials Project energy reference used by Matbench Discovery. Structures and trajectories in ELEMENTA overlapping with the WBM test set were removed during post-processing.

Paper: [Prophet](https://www.kairosmaterials.com/papers/Prophet.pdf)
Code: [kairosmaterial/prophet](https://github.com/kairosmaterial/prophet)

The [submitted discovery predictions](https://doi.org/10.6084/m9.figshare.33691693) use the column name `e_form_per_atom_mathacker-06lx-e150`. The benchmark archive renames it to `e_form_per_atom`, preserves every submitted value, and restores `wbm-5-11353` and `wbm-5-13278` from their supplied relaxed structures and energies using MP2020 corrections. All 256,963 predictions are retained in the archive; the benchmark's central outlier filter excludes those two entries from scoring, leaving the evaluated predictions unchanged.

The submitted runs used Python 3.11 and PyTorch 2.7.1. The registered environment uses Python 3.14 to run the current benchmark code, resolves Prophet's dependencies from its upstream package, and installs the optional CUDA kernels on Linux.
