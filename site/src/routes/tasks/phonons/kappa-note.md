<script lang="ts">
  let { warning = true, details = true } = $props()
</script>

{#if warning}

> ⚠️ To make the κ<sub>SRME</sub> metric more challenging and less susceptible to overfitting, we are working on extending the test set for thermal conductivity prediction from the current size of 103 to at least an order of magnitude more chemically and structurally diverse materials.
> Because it tests the 2nd and 3rd order derivatives of the potential energy surface (PES) and higher derivatives expose even subtle discontinuities in the PES, we believe this to be a stricter and more robust metric for measuring both the utility of ML force fields and the physical accuracy of the PES encoded by an MLFF. It is also interesting because thermal conductivity measurements offer a direct path for benchmarking future MLFFs against experimental data, the highest quality ground truth.
> We invite feedback on this new metric via [GitHub Discussions](https://github.com/janosh/matbench-discovery/discussions/193).

{/if}

{#if details}

> For details on the κ<sub>SRME</sub> modeling task and evaluation method, refer to [arXiv:2408.00755](https://arxiv.org/abs/2408.00755).
> The only difference between the procedure presented by [Póta](https://scholar.google.com/citations?user=5twZvfAAAAAJ), [Ahlawat](https://scholar.google.com/citations?user=-cSYJWUAAAAJ), [Csányi](https://scholar.google.com/citations?user=q39javYAAAAJ), and [Simoncelli](https://tcm.phy.cam.ac.uk/profiles/ms2855), and the results shown here is the relaxation protocol has been simplified and unified for all models (just a single simultaneous cell and site relaxation). See [`matbench_discovery/phonons/thermal_conductivity.py`](https://github.com/janosh/matbench-discovery/blob/7d37186aea19f61806dcc66084a7aaec5ecfbfe0/matbench_discovery/phonons/thermal_conductivity.py) for code to predict 2nd and 3rd order force constants and [`matbench_discovery/metrics/phonons.py`](https://github.com/janosh/matbench-discovery/blob/7d37186aea19f61806dcc66084a7aaec5ecfbfe0/matbench_discovery/metrics/phonons.py) for code to compute κ<sub>SRME</sub>.

{/if}
