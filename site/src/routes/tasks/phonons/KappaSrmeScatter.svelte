<script lang="ts">
  // Per-material κ_SRME vs DFT conductivity for one model, grouped by crystal system.
  import type kappa_data from '$figs/kappa-103-analysis.json.gz'
  import { format_num, sanitize_compact_formula } from 'matterviz'
  import { ScatterPlot } from 'matterviz/plot'
  import type { DataSeries } from 'matterviz/plot'
  import { CRYSTAL_SYSTEM_COLORS, spacegroup_num_to_crystal_sys } from 'matterviz/symmetry'
  import type { CrystalSystem } from 'matterviz/symmetry'
  import { SvelteMap } from 'svelte/reactivity'
  import type { HTMLAttributes } from 'svelte/elements'

  type KappaModelEntry = (typeof kappa_data)[`models`][number]

  interface SrmePoint extends Record<string, unknown> {
    material_id: string
    formula: string
    spg_num: number
    kappa_dft: number
    kappa_ml: number | null
    srme: number
    failures: string
  }

  let { entry, base, ...rest }: HTMLAttributes<HTMLDivElement> & {
    entry: KappaModelEntry
    base: Pick<
      typeof kappa_data,
      `material_ids` | `formulas` | `spg_nums` | `kappa_dft`
    >
  } = $props()

  const failure_labels = [
    [`imag_modes`, `imaginary modes`],
    [`broken_sym`, `broken symmetry`],
    [`max_steps`, `hit max steps`],
  ] as const

  // One series per crystal system lets the legend filter by symmetry.
  let series = $derived.by((): DataSeries<SrmePoint>[] => {
    const by_system = new SvelteMap<CrystalSystem, SrmePoint[]>()
    for (const [idx, mat_id] of base.material_ids.entries()) {
      const [kappa_dft, srme] = [base.kappa_dft[idx], entry.srme[idx]]
      const system = spacegroup_num_to_crystal_sys(base.spg_nums[idx])
      if (kappa_dft === null || srme === null || !system) continue
      const flagged_failures = failure_labels
        .filter(([key]) => entry[key][idx] === true)
        .map(([, label]) => label)
        .join(`, `)
      // Unflagged SRME=2 means the κ calculation crashed.
      const failures =
        flagged_failures || (srme === 2 ? `κ calculation failed` : ``)
      const point: SrmePoint = {
        material_id: mat_id,
        formula: base.formulas[idx],
        spg_num: base.spg_nums[idx],
        kappa_dft,
        kappa_ml: entry.kappa_ml[idx],
        srme,
        failures,
      }
      const points = by_system.get(system) ?? []
      points.push(point)
      by_system.set(system, points)
    }
    return [...by_system.entries()].map(([system, points]) => ({
      x: points.map((point) => point.kappa_dft),
      y: points.map((point) => point.srme),
      metadata: points,
      label: system,
      markers: `points` as const,
      // Censored SRME=2 values render hollow.
      point_style: points.map((point) => ({
        fill: point.srme === 2 ? `transparent` : CRYSTAL_SYSTEM_COLORS[system],
        radius: 5,
        stroke: point.srme === 2 ? CRYSTAL_SYSTEM_COLORS[system] : `white`,
        stroke_width: point.srme === 2 ? 1.5 : 0.5,
      })),
    }))
  })
</script>

<ScatterPlot
  {series}
  x_axis={{ label: `PBE κ<sub>L</sub> (W/mK)`, scale_type: `log`, format: `~s` }}
  y_axis={{ label: `κ<sub>SRME</sub>`, range: [0, 2.05], format: `.1f` }}
  {...rest}
>
  {#snippet tooltip({ metadata })}
    {#if metadata}
      {@const point = metadata as unknown as SrmePoint}
      <strong>{point.material_id}</strong>
      {@html sanitize_compact_formula(point.formula)} (SG {point.spg_num})<br>
      PBE κ: {format_num(point.kappa_dft, `.3~`)} <small>W/mK</small><br>
      {#if point.kappa_ml !== null}
        {entry.label} κ: {format_num(point.kappa_ml, `.3~`)} <small>W/mK</small><br>
      {/if}
      κ<sub>SRME</sub>: {format_num(point.srme, `.3~`)}
      {#if point.failures}<br><em>{point.failures}</em>{/if}
    {/if}
  {/snippet}
</ScatterPlot>
