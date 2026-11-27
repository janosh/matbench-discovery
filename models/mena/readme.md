# Mena

Mena is a family of compact equivariant machine-learning interatomic potentials trained jointly across periodic inorganic crystals, molecular crystals, catalytic surfaces, and open molecules.

- Code: [iatecoding/mena-lab](https://github.com/iatecoding/mena-lab) (`Apache-2.0`)
- Weights: [iatecoding/Mena-LAM-6](https://huggingface.co/iatecoding/Mena-LAM-6/resolve/main/mena-lam6-v1.pt) (`CC-BY-4.0`)
- Predictions archive: [Zenodo 10.5281/zenodo.24461178](https://doi.org/10.5281/zenodo.24461178) & [Hugging Face mirror](https://huggingface.co/datasets/iatecoding/mena-lam-6-matbench-discovery) (`CC-BY-4.0`)
- Paper: [Mena-LAM Technical Report (PDF)](https://github.com/iatecoding/mena-lab/blob/main/docs/Mena-LAM-Technical-Report.pdf)

## Mena-LAM-6

`Mena-LAM-6` (`lam6_v1`, `mena-lam6-v1.pt`, `8,479,848` parameters) is a 3-stage equivariant message-passing potential (`r_cut = 6.0 Å`, `C^3` taper from `5.0` to `6.0 Å`, `env_mod_1152` FFN, continuous coordinate/cloud overlap heads) trained in `48.5` net GPU-hours (`747,692` steps, `Muon + AdamW` with WSD schedule and crystal secant pair loss `P_cry`, strict IEEE FP32 cooldown) on a single NVIDIA RTX PRO 6000 Blackwell GPU across `MPtrj`, `OMat24`, `sAlex`, `OC20`, `MatPES` (PBE & r²SCAN), `OC25`, `OMol25`, `OMC25`, `MAD-1.6`, and `SPICE 2`.

For Materials Project / Matbench Discovery benchmarks, `MenaCalculator` runs in strict IEEE `float32` (`allow_tf32=False`) with `lot="mp"` (factorized `vasp_pbe52` PAW/XC conditioning with composition-based MP2020 `+U` routing: O or F together with V, Cr, Mn, Fe, Co, Ni, Mo, or W) and `add_dispersion=False`.
