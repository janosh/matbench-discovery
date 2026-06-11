<script lang="ts">
  // Failure-mode summary for the kappa-103 task.
  import kappa_data from '$figs/kappa-103-analysis.json.gz'
  import { MODELS } from '$lib'
  import { get_nested_number } from '$lib/metrics'
  import { HeatmapTable } from 'matterviz'
  import type { Label, RowData } from 'matterviz'

  const n_materials = kappa_data.material_ids.length

  const columns: Label[] = [
    { label: `Model`, sticky: true },
    {
      label: `κ<sub>SRME</sub>`,
      description: `Symmetric relative mean error of predicted thermal conductivity (0=perfect, 2=complete failure)`,
      better: `lower`,
      format: `.3~f`,
    },
    {
      label: `κ failed`,
      description:
        `Fraction of the ${n_materials} test materials where the κ calculation failed outright (SRME censored to its max of 2): imaginary phonon modes, broken symmetry during relaxation, or a crashed/NaN κ calculation`,
      better: `lower`,
      format: `.1~%`,
    },
    {
      label: `Imag. modes`,
      description:
        `Fraction of materials with imaginary phonon modes after ML relaxation (unstable predicted structure) - the subset of κ failures with a known physical cause`,
      better: `lower`,
      format: `.1~%`,
    },
    {
      label: `Spectrum W1`,
      description:
        `Mean Wasserstein-1 distance (THz) between ML and DFT phonon frequency spectra. Unlike κ_SRME, this is robust to error compounding in the thermal conductivity calculation and is computed for all materials with phonon frequencies, even those that failed the κ calculation. Values below ~0.02 THz are at the mesh-comparison noise floor (ML spectra are stored on a BZ grid with duplicate zone-boundary points).`,
      better: `lower`,
      format: `.3~f`,
    },
  ]

  const table_data: RowData[] = kappa_data.models.map((entry) => {
    const model = MODELS.find((mod) => mod.model_key === entry.key)
    const srme = model
      ? get_nested_number(model, `metrics.phonons.kappa_103.κ_SRME`)
      : null
    return {
      Model: `<a href="/models/${entry.key}">${entry.label}</a>`,
      'κ<sub>SRME</sub>': srme,
      'κ failed': entry.srme.filter((srme_val) => srme_val === 2).length / n_materials,
      'Imag. modes': entry.imag_modes.filter((flag) => flag === true).length /
        n_materials,
      'Spectrum W1': entry.freq_w1_mean,
    }
  })
</script>

<HeatmapTable
  data={table_data}
  {columns}
  initial_sort={{ column: `κ<sub>SRME</sub>`, direction: `asc` }}
  sort_hint=""
/>
