# Kairos Materials

Prophet is an equivariant message-passing interatomic potential by Kairos Materials
(62.3M parameters, 10 layers, lmax 4) built on Multi-Cutoff Spectral Decomposition,
which models atomic interactions across nested cutoff scales (3/7 A) and recombines
them linearly. It is pretrained on MPtrj, OMat24, sAlex and ELEMENTA, and fine-tuned
for lattice thermal conductivity on Gaussian-displacement fc2/fc3 frames with full-set
replay.

Paper: [Prophet](https://www.kairosmaterials.com/papers/Prophet.pdf)
Code: [kairosmaterial/prophet](https://github.com/kairosmaterial/prophet)
