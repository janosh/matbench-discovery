<script lang="ts">
  import { HeatmapTable, type Column } from 'matterviz/table'
  import { sanitize_compact_formula } from 'matterviz/sanitize'
  import { spacegroup_html } from '#lib/materials.js'
  import type { PageData } from './$types'

  let { data }: { data: PageData } = $props()

  const model_link = (model: { model_key: string; label: string } | null): string =>
    model ? `<a href="/models/${model.model_key}">${model.label}</a>` : `n/a`

  let rows = $derived(
    data.materials.map((material) => ({
      Material: `<a href="/materials/phonons/${material.material_id}">${material.material_id}</a>`,
      Formula: sanitize_compact_formula(material.formula),
      Spacegroup: spacegroup_html(material.spg_num),
      kappa_dft: material.kappa_dft,
      median_srme: material.median_srme,
      median_kappa_ratio: material.median_kappa_ratio,
      n_censored: material.n_censored,
      Best: model_link(material.best),
      Worst: model_link(material.worst),
    })),
  )
  const columns: Column<(typeof rows)[number]>[] = [
    { id: `Material`, label: `Material`, pin: `left` },
    { id: `Formula`, label: `Formula` },
    { id: `Spacegroup`, label: `Spacegroup` },
    {
      id: `kappa_dft`,
      label: `PBE κ<sub>L</sub>`,
      format: `.4~`,
      description: `PhononDB PBE lattice thermal conductivity at 300 K (W/mK)`,
    },
    {
      id: `median_srme`,
      label: `Median κ<sub>SRME</sub>`,
      format: `.3f`,
      better: `lower`,
      description: `Median per-material κ_SRME across models (0 = perfect, 2 = maximum error)`,
    },
    {
      id: `median_kappa_ratio`,
      label: `Median κ<sub>ML</sub>/κ<sub>DFT</sub>`,
      format: `.3~`,
      description: `Median ratio of predicted to PBE κ across models; 1 = exact, below 1 = underestimated`,
    },
    {
      id: `n_censored`,
      label: `N<sub>failed</sub>`,
      format: `d`,
      better: `lower`,
      description: `Models whose κ_SRME was censored to 2 (imaginary modes, broken symmetry or failed κ)`,
    },
    { id: `Best`, label: `Best model`, description: `Lowest κ_SRME on this material` },
    { id: `Worst`, label: `Worst model`, description: `Highest κ_SRME on this material` },
  ]
</script>

<svelte:head>
  <title>Phonon materials • Matbench Discovery</title>
</svelte:head>

<h1>Phonon materials</h1>
<p style="text-align: center; max-width: 50em; margin-inline: auto">
  Per-material results for the {data.materials.length} PhononDB PBE crystals in the
  <a href="/benchmarks/phonons">κ-103 phonon task</a>, across {data.n_models} models. Sort by
  median κ<sub>SRME</sub> to find the materials models struggle with most, then open one to
  compare every model's conductivity and phonon spectrum against DFT.
</p>

<HeatmapTable
  data={rows}
  {columns}
  sort={[{ column: `median_srme`, dir: `desc` }]}
  default_num_format=".3~f"
/>
