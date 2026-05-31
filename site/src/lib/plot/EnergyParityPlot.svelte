<script lang="ts">
  import {
    build_energy_parity_series,
    energy_parity_stats,
    get_energy_parity_point,
    has_energy_parity_model,
    load_energy_parity_base,
    load_energy_parity_model,
    load_wbm_structure,
    structure_popup_placement,
  } from '$lib/energy-parity'
  import type { EnergyKind, EnergyParityBase, EnergyParityModel, EnergyParityPoint, StructurePopupPlacement } from '$lib/energy-parity'
  import type { LoadStatus } from '$lib/asset-loader'
  import type { ModelData } from '$lib/types'
  import { compact_formula, format_num, sanitize_compact_formula } from 'matterviz'
  import { StructurePopup } from 'matterviz/convex-hull'
  import { Spinner } from 'matterviz/feedback'
  import type { AnyStructure } from 'matterviz/structure'
  import { BinnedScatterPlot } from 'matterviz/plot'
  import type { BinnedPointDataFn, BinnedPointPayload, DensePointSeries } from 'matterviz/plot'
  import { onMount, tick, untrack } from 'svelte'
  import type { HTMLAttributes } from 'svelte/elements'

  let {
    model,
    energy_kind,
    title = energy_kind === `e-form` ? `Formation Energies` : `Convex Hull Distance`,
    base_data,
    model_data,
    ...rest
  }: HTMLAttributes<HTMLElement> & {
    model: ModelData
    energy_kind: EnergyKind
    title?: string
    base_data?: EnergyParityBase
    model_data?: EnergyParityModel
  } = $props()

  type EnergyParityPointData = {
    material_id: string
    formula: string
    n_sites: number
    measure_text: string
  }
  const structure_popup_size = { outer_width: 500, view_width: 460, view_height: 340 }
  // small gap lets the popup sit over the (data-free) axis-label padding, so it moves
  // into the side gutter as soon as it clears the data area instead of needing a full
  // popup-width of empty space beside the plot
  const structure_popup_gap = 16
  const loading_spinner_style =
    `--spinner-size: 0.9em; --spinner-border-width: 2px; --spinner-margin: 0`

  let status = $state<LoadStatus>(`idle`)
  let error_message = $state(``)
  let base = $state<EnergyParityBase | undefined>()
  let parity_model = $state<EnergyParityModel | undefined>()
  let selected_point = $state<EnergyParityPoint | null>(null)
  let selected_structure = $state<AnyStructure | null>(null)
  let structure_error = $state(``)
  let structure_loading = $state(false)
  let plot_wrap = $state<HTMLDivElement>()
  let load_id = 0
  let popup_placement = $state<StructurePopupPlacement>({
    side: `left`,
    left: 0,
    top: 0,
  })

  let parity = $derived(
    base && parity_model ? build_energy_parity_series(base, parity_model, energy_kind) : null,
  )
  let stats = $derived(parity ? energy_parity_stats(parity) : null)
  let series = $derived<DensePointSeries[]>(parity
    ? [{
        x: parity.x,
        y: parity.y,
        point_ids: parity.point_ids,
        size_values: parity.size_values,
        label: parity_model?.model_label ?? model.model_name,
        color: model.color ?? `#4dabf7`,
      }]
    : [])
  let axis_label = $derived(
    energy_kind === `e-form` ? `E<sub>form</sub>` : `E<sub>hull dist</sub>`,
  )
  let energy_label = $derived(energy_kind === `e-form` ? `formation energy` : `convex hull distance`)
  let parity_loading_text = $derived(`Loading ${energy_label} parity data...`)
  let structure_loading_text = $derived(`Loading structure for ${energy_label} point...`)
  let x_label = $derived(`PBE ${axis_label}`)
  let y_label = $derived(`${parity_model?.model_label ?? model.model_name} ${axis_label}`)
  // manual min/max loop (not Math.min(...arr)) because WBM x/y arrays are too large
  // to spread as function args without overflowing the call stack
  let extent = $derived.by((): [number, number] => {
    if (!parity) return [-1, 1]
    let [min, max] = [Infinity, -Infinity]
    for (const values of [parity.x, parity.y]) {
      for (const value of values) {
        if (!Number.isFinite(value)) continue
        if (value < min) min = value
        if (value > max) max = value
      }
    }
    if (!Number.isFinite(min) || !Number.isFinite(max)) return [-1, 1]
    if (min === max) return [min - 0.5, max + 0.5]
    const padding = (max - min) * 0.04
    return [min - padding, max + padding]
  })

  async function load_plot_data(
    model_key = model.model_key,
    injected_base = base_data,
    injected_model = model_data,
  ) {
    const current_load_id = ++load_id
    if (!model_key) {
      status = `error`
      error_message = `${model.model_name} has no model key`
      return
    }

    const matching_model = has_energy_parity_model(injected_model, model_key)
      ? injected_model
      : undefined
    if (!has_energy_parity_model(parity_model, model_key)) clear_selection()
    if (injected_base) base = injected_base
    if (matching_model) parity_model = matching_model

    if (base && has_energy_parity_model(parity_model, model_key)) {
      status = `ready`
      return
    }

    status = `loading`
    error_message = ``
    clear_selection()
    try {
      const [base_asset, model_asset] = await Promise.all([
        injected_base ?? base ?? load_energy_parity_base(),
        matching_model ?? load_energy_parity_model(model_key),
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
    const injected_base = base_data
    const injected_model = model_data
    untrack(() => void load_plot_data(model_key, injected_base, injected_model))
  })

  function update_popup_placement() {
    if (!plot_wrap || !selected_point) return
    const rect = plot_wrap.getBoundingClientRect()
    popup_placement = structure_popup_placement({
      viewport_width: globalThis.innerWidth,
      plot_left: rect.left,
      plot_width: rect.width,
      plot_height: rect.height,
      popup_width: structure_popup_size.outer_width,
      gap: structure_popup_gap,
    })
  }

  function clear_selection() {
    selected_point = null
    selected_structure = null
    structure_error = ``
    structure_loading = false
  }

  function energy_parity_point_data(
    point_id: string | number | undefined,
  ): EnergyParityPointData | null {
    const row_idx = Number(point_id)
    if (!base || !Number.isInteger(row_idx)) return null

    const material_id = base.material_ids[row_idx] ?? `unknown-${row_idx}`
    const formula = base.formulas[row_idx] ?? ``
    const n_sites = base.n_sites[row_idx]
    if (typeof n_sites !== `number` || !Number.isFinite(n_sites)) return null
    const measure_text = formula ? `${material_id}\n${compact_formula(formula)}` : material_id
    return { material_id, formula, n_sites, measure_text }
  }

  const point_data: BinnedPointDataFn<Record<string, unknown>, EnergyParityPointData> =
    ({ point }) => energy_parity_point_data(point.point_id)

  function is_energy_parity_point_data(
    point_data_value: unknown,
  ): point_data_value is EnergyParityPointData {
    if (!point_data_value || typeof point_data_value !== `object`) return false
    const { material_id, formula, n_sites, measure_text } =
      point_data_value as Record<keyof EnergyParityPointData, unknown>
    return typeof material_id === `string` &&
      typeof formula === `string` &&
      typeof n_sites === `number` &&
      Number.isFinite(n_sites) &&
      typeof measure_text === `string`
  }

  const point_label_measure_text = ({
    point,
    point_data: point_data_value,
  }: BinnedPointPayload<Record<string, unknown>>): string =>
    is_energy_parity_point_data(point_data_value)
      ? point_data_value.measure_text
      : String(point.point_id ?? ``)

  async function show_structure(point_idx: number) {
    if (!base || !parity_model || !Number.isInteger(point_idx)) return
    const point = get_energy_parity_point(base, parity_model, point_idx, energy_kind)
    if (!point) return
    const selected_row_idx = point.row_idx
    selected_point = point
    selected_structure = null
    structure_error = ``
    structure_loading = true
    void tick().then(update_popup_placement)
    try {
      const structure = await load_wbm_structure(base, point.row_idx, point.material_id)
      if (selected_point?.row_idx === selected_row_idx) selected_structure = structure
    } catch (error) {
      if (selected_point?.row_idx === selected_row_idx) {
        structure_error = error instanceof Error ? error.message : String(error)
      }
    } finally {
      if (selected_point?.row_idx === selected_row_idx) structure_loading = false
    }
  }

  onMount(() => {
    globalThis.addEventListener(`resize`, update_popup_placement)
    return () => {
      globalThis.removeEventListener(`resize`, update_popup_placement)
    }
  })
</script>

<section class="energy-parity-plot" {...rest}>
  <h2 class="toc-exclude">ML vs DFT {title}</h2>

  {#if status === `error`}
    <p class="plot-state" role="alert" style="min-height: 0; margin: 0">{error_message}</p>
  {:else if status !== `ready` || !parity}
    <div class="plot-state">
      <Spinner
        text={parity_loading_text}
        style={loading_spinner_style}
      />
    </div>
  {:else}
    <div class="plot-wrap" bind:this={plot_wrap}>
      {#snippet energy_point_label({ point_data }: BinnedPointPayload<Record<string, unknown>>)}
        {@const label = is_energy_parity_point_data(point_data) ? point_data : null}
        {#if label}
          <span class="point-label-id">{label.material_id}</span>
          {#if label.formula}
            <br><span class="point-label-formula">{@html sanitize_compact_formula(label.formula)}</span>
          {/if}
        {/if}
      {/snippet}

      <BinnedScatterPlot
        {series}
        style="height: 520px"
        x_axis={{ label: x_label, format: `.2f`, range: extent }}
        y_axis={{ label: y_label, format: `.2f`, range: extent }}
        density={{
          color_scale: { type: `log`, scheme: `interpolateMagma` },
          color_bar: { title: `Density` },
        }}
        size_scale={{ radius_range: [2, 18], pick_radius: `auto` }}
        overlays={{
          ref_lines: [{
            x1: extent[0],
            y1: extent[0],
            x2: extent[1],
            y2: extent[1],
            color: `var(--text-color, currentColor)`,
          }],
        }}
        on_point_click={({ point }) => void show_structure(Number(point.point_id))}
        on_density_zoom={clear_selection}
        selected_point_id={selected_point?.row_idx ?? null}
        {point_data}
        point_labels={{
          render: energy_point_label,
          measure_text: point_label_measure_text,
        }}
      >
        {#snippet tooltip({ x, y, x_formatted, y_formatted, point_data: point_data_value })}
          {@const label = is_energy_parity_point_data(point_data_value) ? point_data_value : null}
          {@html x_label}: {x_formatted} <small>eV/atom</small><br>
          {@html y_label}: {y_formatted} <small>eV/atom</small><br>
          MLFF - DFT error: {format_num(y - x, `+.3~`)} <small>eV/atom</small>
          {#if label}
            <br>Points sized by N<sub>atoms</sub>: {format_num(label.n_sites, `.0f`)}
          {/if}
        {/snippet}
      </BinnedScatterPlot>

      {#if stats && Number.isFinite(stats.mae)}
        <div class="plot-annotation">
          MAE = {format_num(stats.mae * 1000, `.3~`)} <small>meV/atom</small><br>
          R<sup>2</sup> = {format_num(stats.r2, `.3~`)}
        </div>
      {/if}

      {#if selected_point}
        {@const point = selected_point}
        <div
          class="popup-anchor {popup_placement.side}"
          style:left={`${popup_placement.left}px`}
          style:top={`${popup_placement.top}px`}
          style:--structure-popup-gap={`${structure_popup_gap}px`}
        >
          {#if selected_structure}
            <StructurePopup
              structure={selected_structure}
              place_right={popup_placement.side === `right`}
              width={structure_popup_size.view_width}
              height={structure_popup_size.view_height}
              stats={{
                formula: point.formula,
              }}
              onclose={clear_selection}
            >
              {#snippet top_left({ formula_html })}
                <strong>{point.material_id}</strong>
                {#if formula_html}
                  ({@html formula_html})<br>
                {:else}
                  <br>
                {/if}
                PBE {@html axis_label}: {format_num(point.x, `.3~`)} <small>eV/atom</small><br>
                {parity_model?.model_label ?? model.model_name} {@html axis_label}:
                {format_num(point.y, `.3~`)} <small>eV/atom</small><br>
                MLFF - DFT error: {format_num(point.y - point.x, `+.3~`)} <small>eV/atom</small>
              {/snippet}
            </StructurePopup>
          {:else}
            <div
              class="structure-status {popup_placement.side}"
              role={structure_error ? `alert` : undefined}
            >
              {#if structure_loading}
                <Spinner
                  text={structure_loading_text}
                  style={loading_spinner_style}
                />
              {:else}
                {structure_error}
              {/if}
            </div>
          {/if}
        </div>
      {/if}
    </div>
  {/if}
</section>

<style>
  .energy-parity-plot {
    margin-block: 2em;
  }
  h2 {
    margin: 1em auto 0.5em;
    text-align: center;
  }
  .plot-state {
    min-height: 180px;
    align-content: center;
    color: var(--muted-text-color, color-mix(in srgb, currentColor 70%, transparent));
    text-align: center;
  }
  .plot-wrap {
    position: relative;
  }
  /* matterviz only sets the tooltip background inline via bg_color, which the binned
     plot omits for per-point tooltips -> give them a readable theme-aware fallback */
  .energy-parity-plot :global(.plot-tooltip) {
    background: var(--tooltip-bg, light-dark(#f5f5f7, #2a2a2e));
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
  }
  .plot-annotation {
    position: absolute;
    right: 2.5em;
    bottom: 4.5em;
    text-align: right;
    font-size: 0.8em;
    line-height: 1.4;
    pointer-events: none;
    background: color-mix(in srgb, var(--surface-bg, rgba(255, 255, 255, 0.6)) 60%, transparent);
    border-radius: 4px;
    padding: 0.1em 0.4em;
  }
  .popup-anchor {
    height: 0;
    position: absolute;
    width: 0;
    z-index: 3;
  }
  .structure-status {
    background: var(--surface-bg, rgba(255, 255, 255, 0.95));
    border: 1px solid var(--border-color, color-mix(in srgb, currentColor 20%, transparent));
    border-radius: 4px;
    box-shadow: 0 16px 24px rgba(0, 0, 0, 0.15);
    min-width: 220px;
    padding: 0.75em 1em;
    position: absolute;
    top: 50%;
    transform: translateY(-50%);
  }
  .structure-status.left {
    right: calc(100% + var(--structure-popup-gap, 1em));
  }
  .structure-status.right {
    left: calc(100% + var(--structure-popup-gap, 1em));
  }
  :global(.energy-parity-plot .structure-popup:not(:hover) .control-buttons) {
    opacity: 0 !important;
    pointer-events: none !important;
  }
  :global(.energy-parity-plot .structure-popup:hover .control-buttons) {
    opacity: 1;
    pointer-events: auto;
  }
  .point-label-id {
    font-weight: 650;
  }
  .point-label-formula {
    font-weight: 400;
  }
</style>
