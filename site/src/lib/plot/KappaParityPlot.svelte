<script lang="ts">
  import {
    as_phonon_dos,
    build_kappa_parity_series,
    has_kappa_parity_model,
    kappa_structure,
    load_kappa_parity_base,
    load_kappa_parity_model,
  } from '$lib/parity/kappa-parity'
  import type {
    KappaParityBase,
    KappaParityModel,
    KappaParityPoint,
    PhononDos,
  } from '$lib/parity/kappa-parity'
  import type { LoadStatus } from '$lib/asset-loader'
  import { parity_diagonal } from '$lib/fig-helpers'
  import { get_nested_number } from '$lib/metrics'
  import type { ModelData } from '$lib/types'
  import { Dos, format_num, Icon, sanitize_compact_formula } from 'matterviz'
  import { Spinner } from 'matterviz/feedback'
  import { ScatterPlot } from 'matterviz/plot'
  import type { DataSeries, RefLine } from 'matterviz/plot'
  import { type CrystalSystem, spacegroup_num_to_crystal_sys } from 'matterviz/symmetry'
  import { untrack } from 'svelte'
  import type { HTMLAttributes } from 'svelte/elements'

  let { model, ...rest }: HTMLAttributes<HTMLElement> & { model: ModelData } = $props()

  let status = $state<LoadStatus>(`idle`)
  let error_message = $state(``)
  let base = $state<KappaParityBase>()
  let parity_model = $state<KappaParityModel>()
  let selected_idx = $state<number | null>(null)
  let load_id = 0

  let model_label = $derived(parity_model?.model_label ?? model.model_name)
  // precomputed phonon metrics: κ_SRME (mode-resolved symmetric relative mean error)
  // and κ_SRE (symmetric relative error of the scalar lattice thermal conductivity)
  let kappa_srme = $derived(get_nested_number(model, `metrics.phonons.kappa_103.κ_SRME`))
  let kappa_sre = $derived(get_nested_number(model, `metrics.phonons.kappa_103.κ_SRE`))
  let parity = $derived(
    base && parity_model ? build_kappa_parity_series(base, parity_model) : null,
  )
  // crystal system (derived from space group) shown in the tooltip for context
  const crystal_sys = (pt: KappaParityPoint): CrystalSystem | null =>
    pt.spacegroup == null ? null : spacegroup_num_to_crystal_sys(pt.spacegroup)
  // shared, multiplicatively-padded range so both axes match and the y=x
  // diagonal runs exactly corner to corner (padding is a factor on a log scale)
  let extent = $derived.by((): [number, number] => {
    if (!parity?.points.length) return [0.1, 100]
    const values = [...parity.x, ...parity.y]
    return [Math.min(...values) / 1.3, Math.max(...values) * 1.3]
  })
  let series = $derived<DataSeries<KappaParityPoint>[]>(
    parity
      ? [
          {
            x: parity.x,
            y: parity.y,
            metadata: parity.points,
            markers: `points`,
            label: model_label,
            size_values: parity.points.map((pt) => pt.n_sites),
            color_values: parity.points.map((pt) => pt.sre),
            point_style: { radius: 6, stroke: `white`, stroke_width: 0.5 },
          },
        ]
      : [],
  )
  const parity_ref_lines: RefLine[] = [parity_diagonal]

  let selected = $derived(
    parity && selected_idx !== null ? (parity.points[selected_idx] ?? null) : null,
  )
  let selected_structure = $derived(
    base && selected ? kappa_structure(base, selected.material_id) : null,
  )
  // three.js stack (~MBs) loads only when a structure is first selected, keeping it
  // out of every page's initial chunk graph
  let Structure = $state<(typeof import('matterviz/structure'))['Structure']>()
  $effect(() => {
    if (selected_structure && !Structure) {
      void import('matterviz/structure').then((mod) => (Structure = mod.Structure))
    }
  })
  let selected_point_ref = $derived(
    selected_idx === null ? null : { series_idx: 0, point_idx: selected_idx },
  )
  let doses = $derived.by((): Record<string, PhononDos> => {
    if (!base || !parity_model || !selected) return {}
    const out: Record<string, PhononDos> = {}
    const dft = as_phonon_dos(base.dft_dos[selected.material_id])
    const ml = as_phonon_dos(parity_model.ml_dos[selected.material_id])
    if (dft) out[`DFT (PBE)`] = dft
    if (ml) out[model_label] = ml
    return out
  })

  async function load_data(model_key = model.model_key) {
    const current_load_id = ++load_id
    selected_idx = null
    if (!model_key || !has_kappa_parity_model(model_key)) {
      status = `error`
      error_message = `${model.model_name} has no κ parity data`
      return
    }
    status = `loading`
    error_message = ``
    try {
      const [base_asset, model_asset] = await Promise.all([
        base ?? load_kappa_parity_base(),
        load_kappa_parity_model(model_key),
      ])
      if (current_load_id !== load_id) return
      base = base_asset
      parity_model = model_asset
      status = `ready`
    } catch (error) {
      if (current_load_id !== load_id) return
      status = `error`
      error_message = error instanceof Error ? error.message : String(error)
    }
  }

  $effect(() => {
    const model_key = model.model_key
    untrack(() => void load_data(model_key))
  })
</script>

<section class="kappa-parity-plot" {...rest}>
  <h2 class="toc-exclude">ML vs DFT Lattice Thermal Conductivity</h2>

  {#if status === `error`}
    <p class="plot-state" role="alert" style="min-height: 0; margin: 0">
      {error_message}
    </p>
  {:else if status !== `ready` || !parity}
    <div class="plot-state">
      <Spinner text="Loading κ parity data..." />
    </div>
  {:else}
    <ScatterPlot
      {series}
      ref_lines={parity_ref_lines}
      style="height: 480px"
      x_axis={{
        label: `PBE κ<sub>L</sub> (W/mK)`,
        scale_type: `log`,
        format: `~s`,
        range: extent,
      }}
      y_axis={{
        label: `${model_label} κ<sub>L</sub> (W/mK)`,
        scale_type: `log`,
        format: `~s`,
        range: extent,
      }}
      size_scale={{ type: `linear`, radius_range: [4, 7] }}
      color_bar={{
        title: `κ<sub>SRE</sub> (${format_num(parity?.points.length ?? 0, `,`)} points)`,
      }}
      selected_point={selected_point_ref}
      on_point_click={({ point }) => {
        if (point.series_idx === 0) selected_idx = point.point_idx
      }}
    >
      {#snippet tooltip({ metadata })}
        {#if metadata}
          {@const pt = metadata as KappaParityPoint}
          {@const sys = crystal_sys(pt)}
          <strong>{pt.material_id}</strong>
          {@html sanitize_compact_formula(pt.formula)}<br />
          PBE κ: {format_num(pt.kappa_dft, `.3~`)} <small>W/mK</small><br />
          {model_label} κ: {format_num(pt.kappa_ml, `.3~`)} <small>W/mK</small><br />
          κ<sub>SRE</sub>: {format_num(pt.sre, `.3~`)}
          {#if sys}<br />{sys}{pt.spacegroup != null ? ` (SG ${pt.spacegroup})` : ``}{/if}
          {#if pt.n_sites != null}<br />{pt.n_sites} atoms{/if}
        {/if}
      {/snippet}

      {#snippet user_content({ width, height })}
        {#if kappa_srme != null || kappa_sre != null}
          {@const style = 'pointer-events: none; overflow: visible'}
          <foreignObject x="0" y="0" {width} {height} {style}>
            <div class="plot-annotation">
              {#if kappa_srme != null}
                κ<sub>SRME</sub> = {format_num(kappa_srme, `.3~`)}<br />
              {/if}
              {#if kappa_sre != null}κ<sub>SRE</sub> = {format_num(kappa_sre, `.3~`)}{/if}
            </div>
          </foreignObject>
        {/if}
      {/snippet}
    </ScatterPlot>

    <p class="plot-note"><small>marker size &propto; atom count</small></p>

    {#if selected}
      <div class="detail-panel">
        <header>
          <span>
            <strong>{selected.material_id}</strong>
            {#if selected.formula}({@html sanitize_compact_formula(
                selected.formula,
              )}){/if}
            &mdash; PBE κ {format_num(selected.kappa_dft, `.3~`)},
            {model_label} κ {format_num(selected.kappa_ml, `.3~`)} <small>W/mK</small>, κ<sub
              >SRE</sub
            >
            {format_num(selected.sre, `.3~`)}
          </span>
          <button
            type="button"
            aria-label="Close"
            title="Close"
            onclick={() => (selected_idx = null)}
          >
            <Icon icon="Cross" />
          </button>
        </header>
        <div class="detail-body">
          {#if selected_structure && Structure}
            {@const style = 'height: 100%; min-height: 360px'}
            <Structure structure={selected_structure} {style} />
          {/if}
          {#if Object.keys(doses).length}
            <Dos {doses} style="height: 360px" padding={{ t: 20, b: 65, r: 10 }} />
          {:else}
            <p class="plot-state">No phonon DOS available for this material.</p>
          {/if}
        </div>
      </div>
    {/if}
  {/if}
</section>

<style>
  .kappa-parity-plot {
    margin-block: 2em;
  }
  h2 {
    margin: 1em auto 0.5em;
    text-align: center;
  }
  .plot-note {
    margin: 0.4em 0 0;
    text-align: center;
    color: var(--muted-text-color, color-mix(in srgb, currentColor 70%, transparent));
  }
  .detail-panel {
    margin-top: 1em;
    border: 1px solid
      var(--border-color, color-mix(in srgb, currentColor 20%, transparent));
    border-radius: 4px;
    overflow: hidden;
  }
  .detail-panel header {
    display: flex;
    gap: 1em;
    align-items: center;
    justify-content: space-between;
    padding: 0.2em 1em;
    background: var(--nav-bg, color-mix(in srgb, currentColor 6%, transparent));
  }
  .detail-panel header button {
    flex: none;
    display: grid;
    place-items: center;
    width: 1.7em;
    height: 1.7em;
    background: color-mix(in srgb, currentColor 10%, transparent);
    border: none;
    border-radius: 50%;
    cursor: pointer;
    opacity: 0;
    padding: 0;
    pointer-events: none;
    transition:
      opacity 0.15s ease,
      background 0.15s ease;
  }
  .detail-panel:is(:hover, :focus-within) header button {
    opacity: 1;
    pointer-events: auto;
  }
  .detail-panel header button :global(svg) {
    display: block;
    height: 1em;
    width: 1em;
  }
  .detail-panel header button:hover {
    background: color-mix(in srgb, currentColor 18%, transparent);
  }
  .detail-body {
    display: grid;
    grid-template-columns: 1fr 1fr;
    column-gap: 1em;
  }
  @media (max-width: 700px) {
    .detail-body {
      grid-template-columns: 1fr;
    }
  }
</style>
