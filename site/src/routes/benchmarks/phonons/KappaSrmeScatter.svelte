<script lang="ts">
  // Per-material κ_SRME vs DFT conductivity for one model, grouped by crystal system.
  // Clicking a point opens that material's cross-model page.
  import { goto } from '$app/navigation'
  import type kappa_data from '#figs/kappa-103-analysis.jsonl'
  import { format_num } from 'matterviz/labels'
  import { sanitize_compact_formula } from 'matterviz/sanitize'
  import { ScatterPlot, type DataSeries } from 'matterviz/plot'
  import { CRYSTAL_SYSTEM_COLORS, spacegroup_to_crystal_sys } from 'matterviz/symmetry'
  import type { CrystalSystem } from 'matterviz/symmetry'
  import type { HTMLAttributes } from 'svelte/elements'

  type KappaModelEntry = (typeof kappa_data)[`models`][number]

  interface SrmePoint extends Record<string, unknown> {
    material_id: string
    formula: string
    spg_num: number
    kappa_dft: number
    kappa_ml: number | null
    kappa_min: number | null
    srme: number
    censored: boolean
    failures: string
  }

  let {
    entry,
    base,
    ...rest
  }: Omit<HTMLAttributes<HTMLDivElement>, `title` | `children`> & {
    entry: KappaModelEntry
    base: Pick<typeof kappa_data, `material_ids` | `formulas` | `spg_nums` | `kappa_dft`>
  } = $props()

  const failure_labels = [
    [`imag_modes`, `imaginary modes`],
    [`broken_sym`, `broken symmetry`],
    [`max_steps`, `hit max steps`],
    [`below_floor`, `κ below half the diffuson minimum`],
  ] as const
  const series_colors: Record<CrystalSystem, string> = {
    ...CRYSTAL_SYSTEM_COLORS,
    hexagonal: `#ab63fa`,
    cubic: `#ef553b`,
  }

  // One series per crystal system lets the legend filter by symmetry.
  let series = $derived.by((): DataSeries<SrmePoint>[] => {
    const by_system = new Map<CrystalSystem, SrmePoint[]>()
    for (const [idx, mat_id] of base.material_ids.entries()) {
      const [kappa_dft, srme] = [base.kappa_dft[idx], entry.srme[idx]]
      const system = spacegroup_to_crystal_sys(base.spg_nums[idx])
      if (kappa_dft === null || srme === null || !system) continue
      const flagged_failures = failure_labels
        .flatMap(([key, label]) => (entry[key][idx] === true ? [label] : []))
        .join(`, `)
      const censored = entry.srme_censored[idx] === true
      const failures = flagged_failures || (censored ? `κ calculation failed` : ``)
      const points = by_system.get(system) ?? []
      points.push({
        material_id: mat_id,
        formula: base.formulas[idx],
        spg_num: base.spg_nums[idx],
        kappa_dft,
        kappa_ml: entry.kappa_ml[idx],
        kappa_min: entry.kappa_min[idx],
        srme,
        censored,
        failures,
      })
      by_system.set(system, points)
    }
    return [...by_system.entries()].map(([system, points]) => ({
      x: points.map((point) => point.kappa_dft),
      y: points.map((point) => point.srme),
      metadata: points,
      label: system,
      markers: `points` as const,
      // Censored values render hollow; a valid SRME can also equal 2.
      point_style: points.map((point) => ({
        fill: point.censored ? `transparent` : series_colors[system],
        radius: 5,
        stroke: point.censored ? series_colors[system] : `white`,
        stroke_width: point.censored ? 1.5 : 0.5,
      })),
    }))
  })
</script>

<ScatterPlot
  {series}
  x_axis={{ label: `PBE κ<sub>L</sub> (W/mK)`, scale_type: `log`, format: `~s` }}
  y_axis={{ label: `κ<sub>SRME</sub>`, range: [0, 2.05], format: `.1f` }}
  on_point_click={({ point }) => {
    const material_id = point.metadata?.material_id
    if (material_id) void goto(`/materials/phonons/${material_id}`)
  }}
  {...rest}
>
  {#snippet tooltip({ metadata })}
    {#if metadata}
      {@const point = metadata as SrmePoint}
      <strong>{point.material_id}</strong>
      {@html sanitize_compact_formula(point.formula)} (SG {point.spg_num})<br />
      PBE κ: {format_num(point.kappa_dft, `.3~`)} <small>W/mK</small><br />
      {#if point.kappa_ml !== null}
        {entry.label} κ: {format_num(point.kappa_ml, `.3~`)} <small>W/mK</small><br />
        {#if point.kappa_min}
          κ / κ<sub>min</sub>: {format_num(point.kappa_ml / point.kappa_min, `.2~`)}
          <small
            title="Agne et al. 2018 diffuson minimum from {entry.label}'s phonon spectrum"
            >(κ<sub>min</sub> {format_num(point.kappa_min, `.2~`)} W/mK)</small
          ><br />
        {/if}
      {/if}
      κ<sub>SRME</sub>: {format_num(point.srme, `.3~`)}
      {#if point.failures}<br /><em>{point.failures}</em>{/if}
    {/if}
  {/snippet}
</ScatterPlot>
