# MD benchmark metric design

Rationale for the molecular-dynamics (finite-temperature) task metrics: each metric
scores how well an MLIP's NVT trajectory reproduces the ab-initio MD references
(DynaMat v1.0, 17 systems, 293–1500 K) on structural, dynamical, and thermodynamic
observables. The public reference artifact intentionally omits energies and forces so it
cannot be used as a supervised validation set. Energy/force RMSEs remain displayable
private-label diagnostics computed only by maintainer infrastructure. This records what
each metric measures, its virtues and flaws, its current discriminating power across
models, and the open questions.

## Metrics at a glance

| metric                  | measures                      | form                                                            | better |
| ----------------------- | ----------------------------- | --------------------------------------------------------------- | ------ |
| `energy_rmse`           | private energy diagnostic     | RMSE of mean-subtracted per-atom energy fluctuations (meV/atom) | lower  |
| `force_rmse`            | private force diagnostic      | RMSE of forces on reference frames (meV/Å)                      | lower  |
| `rdf_error`             | radial structure              | L1 of `g(r)` deviations ÷ distance to ideal gas, %              | lower  |
| `adf_error`             | bond-angle structure          | W1 of bond-angle dist. ÷ distance to `sin θ` background, %      | lower  |
| `vdos_error`            | vibrational dynamics          | W1 of normalized VACF power spectra ÷ σ_ref (band-clipped), %   | lower  |
| `pressure_mae`          | mean stress bias              | `abs(mean(P_ref) − mean(P_pred))` (GPa)                         | lower  |
| `pressure_wasserstein`  | pressure distribution         | W1 between pressure distributions (GPa)                         | lower  |
| `pressure_error` (E_P)  | pressure distribution overlap | L1 histogram non-overlap, %                                     | lower  |
| `combined_score` (CMDS) | aggregate                     | weighted mean of vdos/adf/pressure + speed subscores, [0,1]     | higher |

`rdf_error` is still computed and shown on model detail pages but hidden from the
leaderboard and excluded from CMDS: across models it has Spearman r = 0.91 with
`vdos_error`, 0.92 with `adf_error` and 0.94 with `force_rmse` (the highest force
correlation of any observable), so it double-counts structural accuracy. Dropping it
from CMDS changes model rankings by at most 1 place (Spearman 0.9994 between the
3- and 4-component scores over N=34 models).

**Paper lineage.** `rdf_error`, the pressure metrics (`pressure_mae`, `pressure_error`) and
`vdos_error` track the reference paper's MD observables. `adf_error` and `combined_score`
(CMDS) are Matbench *additions* not in the paper. `vdos_error` keeps the paper's observable
but swaps its L1 spectral overlap for a W1/σ_ref distance (see below). Single-point
energy/force RMSEs require private supervised AIMD labels and are intentionally excluded
from CMDS.

## Discriminating power

(N = 33 models at the time of this snapshot; correlation numbers elsewhere use the
current N = 34 cohort.)

| metric                 | min  | max   | mean | CV (std/mean) | note                                             |
| ---------------------- | ---- | ----- | ---- | ------------- | ------------------------------------------------ |
| `energy_rmse`          | 0.49 | 6.06  | 1.18 | **0.83**      | private-label diagnostic                         |
| `force_rmse`           | 43.3 | 202.2 | 99.6 | 0.37          | private-label diagnostic                         |
| `rdf_error`            | 3.09 | 9.66  | 4.58 | 0.30          | moderate                                         |
| `adf_error`            | 1.28 | 6.27  | 1.97 | 0.47\*        | \*outlier-driven (one model 6.3%, rest 1.3–2.8%) |
| `vdos_error`           | 21.6 | 35.6  | 24.2 | 0.13          | low; ~2x its 0.06 under L1 overlap               |
| `pressure_mae`         | 0.70 | 1.21  | 0.84 | 0.13          | low                                              |
| `pressure_wasserstein` | 0.71 | 1.21  | 0.85 | 0.13          | low                                              |
| `pressure_error`       | 42.2 | 73.8  | 54.8 | 0.16          | low–moderate                                     |
| `combined_score`       | 0.71 | 0.83  | 0.79 | 0.04          | compressed (it is a [0,1] score)                 |

Cross-metric Pearson `r`: `adf↔rdf` 0.94, `rdf↔vdos` 0.95, `adf↔force` 0.81;
`pressure_error↔force` 0.75, `↔rdf` 0.70, `↔vdos` 0.64. The structural/dynamical metrics
are highly mutually correlated (and with force RMSE) — notably W1-vDOS is **more** RDF-
correlated (0.95) than L1-vDOS was (0.90); pressure is the most independent.

## Distance metric: L1 overlap vs Wasserstein-1 (W1)

Several metrics compare two area-normalized 1-D distributions (spectra/histograms). Two
choices of distance, with very different behavior:

- **L1 (overlap)** = `∫|f_ref − f_pred|` — a *vertical*, pointwise difference; the "%
  non-overlapping area". Only credits mass that lands in the *same* bin.
- **W1 (earth-mover)** = `∫|CDF_ref − CDF_pred|` — a *horizontal* transport distance; the
  work (mass × distance) to reshape one distribution into the other.

**Key behavioral difference — saturation/inversion vs monotonic growth.** Synthetic test
(reference bands at 5 & 15 THz, both shifted by Δ):

| shift Δ (THz) | L1 (%)   | W1 (THz) |
| ------------- | -------- | -------- |
| 0.25          | 14.5     | 0.25     |
| 1.0           | 53.4     | 1.00     |
| 2.0           | 84.9     | 2.00     |
| 5.0           | 99.9     | 5.00     |
| 10.0          | **57.6** | 10.00    |

W1 is exactly monotonic. **L1 saturates by ~5 THz and then *decreases*** — a grossly
shifted distribution scores *better* than a mildly shifted one because the displaced peak
accidentally re-overlaps the reference. So L1 not only loses resolution on large shifts,
it can **invert rankings**. This matters because the dominant structural/dynamical failure
mode (bond-length, bond-angle, phonon softening/hardening) is a *shift*.

W1 virtues/flaws:

- **+** grows with how far mass moved → resolves shifts that L1 cannot; also captures
  broadening monotonically (it is not blind to peak shape).
- **−** unbounded (axis units), so it needs normalization (see below).
- **−** distance-weighting makes it sensitive to far-tail mass: a small spurious peak far
  out, or estimation noise, is penalized heavily and even ranked by *where* the noise
  lands. **Mitigation:** clip the distribution to the physically meaningful band before
  computing W1 (e.g. to 99.9% of reference power); in synthetic tests this drives the
  tail-noise contribution to ~0 while preserving shift sensitivity.

**Normalization of W1 to a bounded, comparable score:**

- where a "featureless" reference exists, divide by the distance to it — the bounded
  "fraction of the way to random" (ADF uses `W1(ref,pred)/W1(ref, sin θ)`; RDF's analog is
  the ideal gas `g = 1`).
- otherwise divide by the reference's own scale, `W1/σ_ref` (its spectral/distribution
  std) — dimensionless and scale-invariant, with a floor for near-degenerate references.

L1 virtues/flaws:

- **+** bounded in [0, 100%] and robust: each unit of misplaced mass contributes the same
  regardless of location (insensitive to far-tail noise).
- **−** saturates and can invert rankings on shifts (above); blind to *how far* mass moved.

## Per-metric notes

- **`energy_rmse`** — private-label mean-subtracted fluctuation RMSE; the mean subtraction
  removes the absolute-energy offset, which is not comparable across systems with different
  DFT energy zeros (all-electron vs PAW). It also makes the metric **functional-agnostic**:
  by scoring energy *fluctuations* rather than absolute energies, it stays comparable across
  a change of reference DFT functional (e.g. the planned PBE → r2SCAN migration), where
  absolute energies shift but the ensemble-relevant fluctuations are preserved. Excluded
  from CMDS because the public reference artifact omits energies.
- **`force_rmse`** — private-label force diagnostic; improves systematically across model
  generations and correlates strongly with the structural/dynamical observables. Excluded
  from CMDS because the public reference artifact omits forces.
- **`rdf_error`** — L1 of `g(r)` deviations normalized by distance to the ideal gas; uses a
  minimum-image radial cutoff. Flaws: global all-pair pooling can hide minority-environment
  failures, and L1 saturates on bond-length shifts. A W1 form (on the pair-distance density
  `∝ r²ρg(r)`, normalized by the ideal-gas distance) is the shift-sensitive alternative but
  requires reformulating away from `g(r)` (which is a ratio, not a density).
- **`adf_error`** — angular complement to the RDF over species-aware covalent-radius bonded
  triplets, using **W1** (deliberately, to resolve angular shifts like tetrahedral 109.5° →
  square-planar 90° without saturating). Flaws: (1) global triplet pooling weights common
  coordination environments by raw count, hiding minority failures; (2) ~0.94 correlated
  with RDF → limited *independent* signal; (3) spread is driven by a single outlier. A
  per-unordered-species-triplet form (support-weighted, missing reference triplet → 100%)
  is a candidate to surface minority-environment errors and decorrelate from RDF.
- **`vdos_error`** — **W1** between Hann-windowed VACF power spectra (area-normalized, on the
  common grid up to the smaller Nyquist), normalized by the reference spectral std `σ_ref`
  and capped at 100%, after clipping the reference tail beyond 99.9% of its cumulative power.
  Chosen over the former L1 overlap because the dominant failure mode is a *shift* where L1
  saturates/inverts (quantified in "W1 vs L1 for vDOS" below). Sensitive to the saved-frame
  interval `dt_fs` (velocities ∝ 1/dt, so a wrong `dt_fs` rescales
  the frequency axis); guarded by an equipartition-temperature sanity check on the reference.
- **pressure** — `P = −tr(σ)/3`. Three views:
  - `pressure_mae = |mean(P_ref) − mean(P_pred)|` (mean-stress bias). Frame-by-frame pairing
    is *not* used: reference and MLIP trajectories are independently thermalized, so paired
    instantaneous pressures decorrelate and a paired MAE measures fluctuation noise, not
    model quality.
  - `pressure_wasserstein` — W1 of the pressure distributions (GPa); robust to frame
    decorrelation, captures both mean offset and width, and does **not** saturate.
  - `pressure_error` (E_P) — L1 histogram non-overlap. Flaw: saturates — all distributions
    with little/no overlap return ≈100%, so it cannot discriminate among the worst models.
  Pressure is the least correlated with the other observables and the hardest to predict
  (often attributed to noisy stress labels in large training sets).

## Combined score (CMDS)

CMDS is a weighted mean of four subscores in [0,1], higher is better: the observables
`vdos_error`, `adf_error` and `pressure_error` (each with subscore `1 − error/100`) plus
**speed** — the summed NVT rollout wall time `run_time_sec`, log-clamped from 9e3 s
(subscore 1) to 3e5 s (subscore 0), bounds frozen at the mid-2026 field (33 models, one
H200 per system), mirroring the CDS speed pillar. All four are **equally weighted at 1/4
by default**, which sits the radar knob dead center and is trivially knob-expressible
(opposite-corner products 1/16 == 1/16, see `DEFAULT_CDS_CONFIG`). The square knob can't
express a vDOS-dominant default at this speed weight — the diagonal-opposite product would
force ADF·pressure too large — so equal weighting is the natural pick; users can reweight
live on the MD task page. Models without a recorded `run_time_sec` get no CMDS unless the
speed weight is zeroed.
`rdf_error` was dropped from the aggregate (and hidden from the leaderboard, kept on model
detail pages): across N=34 models it is Spearman-correlated 0.91 with `vdos_error`, 0.92
with `adf_error` and 0.94 with `force_rmse`, so it double-counted structural accuracy while
adding a column. The 3-observable CMDS agrees with the old 4-observable one at Spearman
0.9994 (max rank shift: 1 place).

**Pressure dominates the aggregate.** Because the components live on very different raw
scales (adf ≈ 2%, vdos ≈ 24%, pressure ≈ 55%), an equal-weight mean is not an
equal-influence mean:

| component        | mean error | fraction of mean CMDS loss |
| ---------------- | ---------- | -------------------------- |
| `adf_error`      | 2.0        | 2%                         |
| `vdos_error`     | 24.2       | 30%                        |
| `pressure_error` | 54.8       | **68%**                    |

Empirically `pressure_error` alone explains **R² ≈ 0.9** of the CMDS ranking — CMDS is
largely a pressure ranking. Aggregation options considered:

- **raw-error mean** (kept): simple; lets the largest-magnitude (hardest) error
  dominate.
- **Lehmer mean** `Σx²/Σx`: amplifies the dominant error further (more pressure-weighted).
- **per-component calibration**: map each error to a comparable subscore (noise floor +
  bad-baseline) before averaging, so components contribute equally.

**Per-component calibration: measured, and rejected for live use.** Leave-one-out
experiment on the N=34 model-level metrics — recompute each calibrated combined score with
one model removed from the cohort and record the worst rank shift among the remaining
models:

| aggregate           | spread (std) | max rank shift on removing 1 model |
| ------------------- | ------------ | ---------------------------------- |
| raw mean (current)  | 0.038        | **0**                              |
| min-max calibrated  | 0.201        | 4                                  |
| z-score calibrated  | 0.150        | 3                                  |

Calibration buys 4–5× more spread but makes every published score a function of the model
roster: adding or removing a single model reshuffles existing models by up to 4 ranks
without any of them changing. The raw mean is perfectly stable under roster changes. If
more spread is ever wanted, the stable option is **normalization scales frozen at the
v1.0 cohort** (fixed min/max per component, documented and never recomputed), which keeps
scores roster-independent by construction.

Decision: keep the **simple weighted mean of raw-error subscores** as a `[0,1]` score (no
roster-dependent calibration). Naming/semantics: it is a *score* (higher = better),
distinct from the component *errors* (lower = better).

CMDS is **not stored** in model YAMLs (only its components are): the site computes it on
the fly like CPS, so a formula change can never leave stale scores behind, and users can
reweight the components live via the radar-chart UI on the MD task page.

## Decisions (implemented)

- **W1 over L1 for shift-dominated distributions** (with band-clipping + a principled
  normalization) — W1 is monotonic where L1 saturates/inverts. Adopt per-metric on
  evidence, not blanket.
- **`vdos_error` = W1/σ_ref (band-clipped)**, replacing the L1 overlap for ~2× the
  discriminating spread and monotonicity (no rank inversion on phonon shifts). NB: unlike
  the experiment's raw-density variant, the shipped (mass-weighted) metric is *not* less
  redundant with RDF — see results below.
- **`pressure_mae` = mean-stress bias**, not a frame-paired MAE (independent trajectories
  decorrelate).
- **CMDS** stays a simple weighted `[0,1]` mean score (observables + log-scaled speed,
  no roster-dependent calibration).
- **`adf_error`** uses W1 + background normalization; minority-environment / per-triplet
  redesign is a candidate, gated on whether it adds signal independent of RDF.
- Input validation (finite/non-negative) on all distribution-error functions.

## Result: W1 vs L1 for vDOS (36 models)

Empirical test recomputing both metrics on identical spectra for every cached rollout:

| variant                    | CV (spread) | Spearman ρ vs L1 | Pearson r vs RDF | vs force |
| -------------------------- | ----------- | ---------------- | ---------------- | -------- |
| L1 overlap (former)        | 0.059       | —                | 0.90             | 0.89     |
| W1 raw (THz)               | 0.048       | 0.73             | —                | —        |
| **W1/σ_ref, band-clipped** | **0.121**   | 0.91             | 0.71\*           | 0.71\*   |

\*raw-density experiment variant. The **shipped** metric adds trapezoidal mass-weighting and
is **0.95** vs RDF / 0.88 vs force on the 33 committed models (discriminating-power table) —
i.e. as redundant as L1; the redundancy reduction below did not carry to production.

- **~2× discriminating spread**: W1/σ_ref CV 0.121 vs L1 0.059. The σ_ref normalization is
  essential — *raw* W1 in THz (CV 0.048) is no better than L1.
- **Redundancy — did *not* survive to production**: the raw-density W1/σ here dropped the RDF
  correlation from L1's ≈0.90 to ≈0.71, but the shipped mass-weighted metric is back to ≈0.95
  (above). So W1's retained edge over L1 is the spread and monotonicity, not independence.
- **Robust**: identical to 4 decimal places under Gaussian smoothing (corr 1.0000), so the
  far-tail-noise concern does not bite in practice once band-clipped.
- **Caveat**: W1 reshuffles the *top* of the ranking (e.g. it promotes some older Tier-2
  models that L1 rates mid-pack; worst models agree). It is not merely "same order, more
  spread," so promoting W1 to the headline vDOS metric should include a physical spot-check
  of the spectra of the promoted/demoted models.

**Force-RMSE cross-check (20 models with both vdos-variant and YAML force metrics).** How
well each vDOS variant's model-level mean tracks force quality:

| variant                    | CV (spread) | Pearson r vs force | Spearman ρ vs force |
| -------------------------- | ----------- | ------------------ | ------------------- |
| L1 overlap                 | 0.059       | 0.88               | 0.88                |
| W1 raw (THz)               | 0.048       | 0.57               | 0.63                |
| **W1/σ_ref, band-clipped** | **0.121**   | 0.83               | **0.93**            |

W1/σ_ref band-clipped has both the highest spread *and* the highest rank agreement with
force RMSE; L1 edges it only on Pearson (linear) correlation. Raw unnormalized W1 is worse
on every axis — the σ_ref normalization does the work. This confirms the choice of
W1/σ_ref (band-clipped) as the shipped vDOS metric. The two rank models nearly identically
overall (L1 vs W1/σ Spearman 0.91); the differences concentrate at the top of the table.

## Open questions / in-progress

- **vDOS top-of-ranking spot-check**: W1/σ_ref is now the headline metric; a physical
  spot-check of the spectra of the models it promotes/demotes vs L1 is still advisable.
- **ADF**: per-species-triplet redesign to surface minority-environment errors (with RDF
  out of CMDS, ADF is the sole structural component, strengthening the case).
- **RDF**: optional W1 (pair-distance density) and species-resolved components, if it is
  ever to earn its way back into CMDS as an independent signal.
- **Energy-distribution W1**: a bounded 0–1 alternative to the energy-fluctuation RMSE
  (W1 between the sampled energy distributions, normalized like the other observables).
  Moot for the public leaderboard while energies are private-label diagnostics, but the
  natural form if energy metrics ever become public (e.g. via a released labeled subset).
