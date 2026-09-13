# Kairos Materials

Prophet is an equivariant message-passing interatomic potential by Kairos Materials (62.3M parameters, 10 layers, lmax 4) built on Multi-Cutoff Spectral Decomposition, which models atomic interactions across nested cutoff scales (3/7 A) and recombines them linearly. It is pretrained on MPtrj, OMat24, sAlex and ELEMENTA, and fine-tuned for lattice thermal conductivity on Gaussian-displacement fc2/fc3 frames with full-set replay.

Paper: [Prophet](https://www.kairosmaterials.com/papers/Prophet.pdf)
Code: [kairosmaterial/prophet](https://github.com/kairosmaterial/prophet)

The [submitted discovery predictions](https://doi.org/10.6084/m9.figshare.33691693) use the column name `e_form_per_atom_mathacker-06lx-e150`. The benchmark archive renames it to `e_form_per_atom`, preserves every submitted value, and restores `wbm-5-11353` and `wbm-5-13278` from their supplied relaxed structures and energies using MP2020 corrections. All 256,963 predictions are retained in the archive; the benchmark's central outlier filter excludes those two entries from scoring, leaving the evaluated predictions unchanged.

The submitted runs used Python 3.11 and PyTorch 2.7.1. The registered environment uses Python 3.14 to run the current benchmark code, resolves Prophet's dependencies from its upstream package, and installs the optional CUDA kernels on Linux.
