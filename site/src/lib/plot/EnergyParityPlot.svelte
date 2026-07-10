<script lang="ts">
  import {
    build_energy_parity_series,
    energy_parity_stats,
    get_energy_parity_point,
    load_energy_parity_base,
    load_energy_parity_model,
    load_wbm_structure,
    structure_popup_placement,
  } from '$lib/parity/energy-parity'
  import type {
    EnergyKind,
    EnergyParityBase,
    EnergyParityModel,
    EnergyParityPoint,
    StructurePopupPlacement,
  } from '$lib/parity/energy-parity'
  import { get_error_message, type LoadStatus } from '$lib/asset-loader'
  import type { ModelData } from '$lib/types'
  import { compact_formula, format_num, sanitize_compact_formula } from 'matterviz'
  import { Spinner } from 'matterviz/feedback'
  import type { AnyStructure } from 'matterviz/structure'
  import { BinnedScatterPlot } from 'matterviz/plot'
  import type {
    BinnedOverlaysConfig,
    BinnedPointDataFn,
    BinnedPointPayload,
    DensePointSeries,
  } from 'matterviz/plot'
  import { onMount, tick, untrack } from 'svelte'
  import type { HTMLAttributes } from 'svelte/elements'

  let {
    model,
    energy_kind,
    onstatus,
    ...rest
  }: HTMLAttributes<HTMLElement> & {
    model: ModelData
    energy_kind: EnergyKind
    // notified on load-status changes, e.g. to show a spinner in a tab bar
    onstatus?: (status: LoadStatus) => void
  } = $props()

  type EnergyParityPointData = {
    material_id: string
    formula: string
    n_sites: number
    measure_text: string
  }
  type EnergyParityMetadata = Record<string, unknown>
  type EnergyParityPayload = BinnedPointPayload<
    EnergyParityMetadata,
    EnergyParityPointData
  >
  const EnergyParityScatter = BinnedScatterPlot<
    EnergyParityMetadata,
    EnergyParityPointData
  >
  const structure_popup_size = { outer_width: 500, view_width: 460, view_height: 340 }
  const parity_overlays: BinnedOverlaysConfig = {
    ref_lines: [
      {
        type: `diagonal`,
        slope: 1,
        intercept: 0,
        style: { color: `var(--text-color, currentColor)` },
      },
    ],
  }
  // small gap lets the popup sit over the (data-free) axis-label padding, so it moves
  // into the side gutter as soon as it clears the data area instead of needing a full
  // popup-width of empty space beside the plot
  const structure_popup_gap = 16
  const loading_spinner_style = `--spinner-size: 0.9em; --spinner-border-width: 2px; --spinner-margin: 0`

  // Wait for loading UI to paint; timeout avoids throttled rAF in background tabs.
  const wait_for_loading_paint = (): Promise<void> =>
    new Promise((resolve) => {
      const timeout_id = globalThis.setTimeout(resolve, 100)
      globalThis.requestAnimationFrame(() =>
        globalThis.requestAnimationFrame(() => {
          globalThis.clearTimeout(timeout_id)
          resolve()
        }),
      )
    })

  let status = $state<LoadStatus>(`idle`)
  $effect(() => onstatus?.(status))
  let error_message = $state(``)
  let base = $state<EnergyParityBase | undefined>()
  let parity_model = $state<EnergyParityModel | undefined>()
  let selected_point = $state<EnergyParityPoint | null>(null)
  let selected_structure = $state<AnyStructure | null>(null)
  let structure_error = $state(``)
  let structure_loading = $state(false)
  let plot_wrap = $state<HTMLElement>()
  let load_id = 0
  // three.js stack (~MBs) loads only when a structure is first clicked, keeping it
  // out of every page's initial chunk graph
  let StructurePopup =
    $state<(typeof import('matterviz/convex-hull'))['StructurePopup']>()
  let popup_placement = $state<StructurePopupPlacement>({ side: `left`, left: 0, top: 0 })

  let model_label = $derived(parity_model?.model_label ?? model.model_name)
  let parity = $derived(
    base && parity_model
      ? build_energy_parity_series(base, parity_model, energy_kind)
      : null,
  )
  let stats = $derived(parity ? energy_parity_stats(parity) : null)
  let series = $derived<DensePointSeries[]>(
    parity
      ? [
          {
            x: parity.x,
            y: parity.y,
            point_ids: parity.point_ids,
            size_values: parity.size_values,
            label: model_label,
            color: model.color ?? `#4dabf7`,
          },
        ]
      : [],
  )
  let axis_label = $derived(
    energy_kind === `e-form` ? `E<sub>form</sub>` : `E<sub>hull dist</sub>`,
  )
  let energy_label = $derived(
    energy_kind === `e-form` ? `formation energy` : `convex hull distance`,
  )
  let x_label = $derived(`PBE ${axis_label}`)
  let y_label = $derived(`${model_label} ${axis_label}`)
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

  async function load_plot_data(model_key: string | undefined) {
    const current_load_id = ++load_id
    if (!model_key) {
      status = `error`
      error_message = `${model.model_name} has no model key`
      return
    }

    if (base && parity_model?.model_key === model_key) {
      status = `ready`
      return
    }

    status = `loading`
    error_message = ``
    clear_selection()
    try {
      await wait_for_loading_paint()
      if (current_load_id !== load_id) return
      const [base_asset, model_asset] = await Promise.all([
        load_energy_parity_base(),
        load_energy_parity_model(model_key),
      ])
      if (current_load_id !== load_id) return
      base = base_asset
      parity_model = model_asset
      // Keep status=loading while the expensive derived series and plot DOM render.
      await tick()
      if (current_load_id !== load_id) return
      status = `ready`
    } catch (error) {
      if (current_load_id !== load_id) return
      status = `error`
      error_message = get_error_message(error)
    }
  }

  $effect(() => {
    const model_key = model.model_key
    untrack(() => void load_plot_data(model_key))
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
    const measure_text = formula
      ? `${material_id}\n${compact_formula(formula)}`
      : material_id
    return { material_id, formula, n_sites, measure_text }
  }

  const point_data: BinnedPointDataFn<EnergyParityMetadata, EnergyParityPointData> = ({
    point,
  }) => energy_parity_point_data(point.point_id)

  async function show_structure(point_idx: number) {
    if (!base || !parity_model || !Number.isInteger(point_idx)) return
    const point = get_energy_parity_point(base, parity_model, point_idx, energy_kind)
    if (!point) return
    const selection_is_current = () => selected_point?.row_idx === point.row_idx
    selected_point = point
    selected_structure = null
    structure_error = ``
    structure_loading = true
    void tick().then(update_popup_placement)
    try {
      const [structure] = await Promise.all([
        load_wbm_structure(base, point.row_idx, point.material_id),
        StructurePopup ??
          import('matterviz/convex-hull').then((mod) => {
            StructurePopup = mod.StructurePopup
          }),
      ])
      if (selection_is_current()) selected_structure = structure
    } catch (error) {
      if (selection_is_current()) {
        structure_error = get_error_message(error)
      }
    } finally {
      if (selection_is_current()) structure_loading = false
    }
  }

  onMount(() => {
    globalThis.addEventListener(`resize`, update_popup_placement)
    return () => globalThis.removeEventListener(`resize`, update_popup_placement)
  })

  // matterviz auto-places the density colorbar in whichever corner least occludes
  // data, so the MAE/R² annotation claims the diagonally opposite corner to
  // guarantee the two never overlap. Insets clear the axes + their tick labels.
  // colorbar_class is owned by this file and injected through matterviz's public
  // color_bar props, so placement never depends on matterviz-internal class names.
  const colorbar_class = `density-color-bar`
  const colorbar_selector = `.${colorbar_class}`
  const annotation_insets = {
    top_left: `2.5em auto auto 7em`,
    top_right: `2.5em 2em auto auto`,
    bottom_left: `auto auto 5em 7em`,
    bottom_right: `auto 2em 5em auto`,
  }
  let annotation_inset = $state(annotation_insets.bottom_right)

  function place_annotation_opposite_colorbar() {
    const bar = plot_wrap?.querySelector(colorbar_selector)?.getBoundingClientRect()
    const wrap = plot_wrap?.getBoundingClientRect()
    // zero-width wrap = plot in a hidden tab; its rects would misplace the annotation
    if (!bar || !wrap?.width) return
    const vert = bar.top + bar.height / 2 < wrap.top + wrap.height / 2 ? `bottom` : `top`
    const horiz = bar.left + bar.width / 2 < wrap.left + wrap.width / 2 ? `right` : `left`
    annotation_inset = annotation_insets[`${vert}_${horiz}`]
  }

  $effect(() => {
    if (status !== `ready` || !plot_wrap) return
    // the colorbar mounts late and moves on zoom/resize (all via inline-style
    // updates), so watch mutations involving it instead of enumerating triggers.
    // Cheap filter keeps tooltip style churn from forcing layout on every mousemove.
    const involves_colorbar = (node: Node) =>
      node instanceof HTMLElement &&
      (node.closest(colorbar_selector) ?? node.querySelector(colorbar_selector)) != null
    const observer = new MutationObserver((mutations) => {
      if (
        mutations.some((mut) => [mut.target, ...mut.addedNodes].some(involves_colorbar))
      ) {
        place_annotation_opposite_colorbar()
      }
    })
    observer.observe(plot_wrap, {
      childList: true,
      subtree: true,
      attributes: true,
      attributeFilter: [`style`],
    })
    place_annotation_opposite_colorbar()
    return () => observer.disconnect()
  })
</script>

<!-- the plot has no visible heading of its own: the model page's tab bar acts as its
title, so label the section for screen readers instead -->
<section
  class="energy-parity-plot"
  aria-label="ML vs DFT {energy_label} parity plot"
  bind:this={plot_wrap}
  {...rest}
>
  {#if status === `error`}
    <p class="plot-state" role="alert" style="min-height: 0; margin: 0">
      {error_message}
    </p>
  {:else if !parity || parity_model?.model_key !== model.model_key}
    <div class="plot-state">
      <Spinner
        text="Loading {energy_label} parity data..."
        style={loading_spinner_style}
      />
    </div>
  {:else}
    {#snippet energy_point_label({ point_data }: EnergyParityPayload)}
      {#if point_data}
        <span class="point-label-id">{point_data.material_id}</span>
        {#if point_data.formula}
          <br /><span class="point-label-formula"
            >{@html sanitize_compact_formula(point_data.formula)}</span
          >
        {/if}
      {/if}
    {/snippet}

    <EnergyParityScatter
      {series}
      style="height: 520px"
      x_axis={{ label: x_label, format: `.2f`, range: extent }}
      y_axis={{ label: y_label, format: `.2f`, range: extent }}
      density={{
        color_scale: { type: `log`, scheme: `interpolateMagma` },
        color_bar: { title: `Density`, class: colorbar_class },
      }}
      size_scale={{ radius_range: [2, 18], pick_radius: `auto` }}
      overlays={parity_overlays}
      on_point_click={({ point }) => void show_structure(Number(point.point_id))}
      on_density_zoom={clear_selection}
      selected_point_id={selected_point?.row_idx ?? null}
      {point_data}
      point_labels={{
        render: energy_point_label,
        measure_text: ({ point, point_data }: EnergyParityPayload) =>
          point_data?.measure_text ?? String(point.point_id ?? ``),
      }}
    >
      {#snippet tooltip({ x, y, x_formatted, y_formatted, point_data })}
        {@html x_label}: {x_formatted} <small>eV/atom</small><br />
        {@html y_label}: {y_formatted} <small>eV/atom</small><br />
        MLFF - DFT error: {format_num(y - x, `+.3~`)} <small>eV/atom</small>
        {#if point_data}
          <br />Points sized by N<sub>atoms</sub>: {format_num(point_data.n_sites, `.0f`)}
        {/if}
      {/snippet}

      {#snippet children()}
        {#if stats && Number.isFinite(stats.mae)}
          <div class="plot-annotation" style:inset={annotation_inset}>
            MAE = {format_num(stats.mae * 1000, `.3~`)} <small>meV/atom</small><br />
            R<sup>2</sup> = {format_num(stats.r2, `.3~`)}
          </div>
        {/if}
      {/snippet}
    </EnergyParityScatter>

    {#if selected_point}
      {@const point = selected_point}
      <div
        class="popup-anchor {popup_placement.side}"
        style:left="{popup_placement.left}px"
        style:top="{popup_placement.top}px"
        style:--structure-popup-gap="{structure_popup_gap}px"
      >
        {#if selected_structure && StructurePopup}
          <StructurePopup
            structure={selected_structure}
            place_right={popup_placement.side === `right`}
            width={structure_popup_size.view_width}
            height={structure_popup_size.view_height}
            stats={{ formula: point.formula }}
            onclose={clear_selection}
          >
            {#snippet top_left({ formula_html })}
              <strong>{point.material_id}</strong>
              {#if formula_html}
                ({@html formula_html})<br />
              {:else}
                <br />
              {/if}
              PBE {@html axis_label}: {format_num(point.x, `.3~`)}
              <small>eV/atom</small><br />
              {parity_model?.model_label ?? model.model_name}
              {@html axis_label}:
              {format_num(point.y, `.3~`)} <small>eV/atom</small><br />
              MLFF - DFT error: {format_num(point.y - point.x, `+.3~`)}
              <small>eV/atom</small>
            {/snippet}
          </StructurePopup>
        {:else}
          <div
            class="structure-status {popup_placement.side}"
            role={structure_error ? `alert` : undefined}
          >
            {#if structure_loading}
              <Spinner
                text="Loading structure for {energy_label} point..."
                style={loading_spinner_style}
              />
            {:else}
              {structure_error}
            {/if}
          </div>
        {/if}
      </div>
    {/if}
  {/if}
</section>

<style>
  .energy-parity-plot {
    position: relative;
  }
  /* matterviz only sets the tooltip background inline via bg_color, which the binned
     plot omits for per-point tooltips -> give them a readable theme-aware fallback */
  .energy-parity-plot :global(.plot-tooltip) {
    background: var(--tooltip-bg, light-dark(#f5f5f7, #2a2a2e));
    box-shadow: 0 4px 12px var(--shadow, rgba(0, 0, 0, 0.15));
  }
  .popup-anchor {
    height: 0;
    position: absolute;
    width: 0;
    z-index: 3;
  }
  .structure-status {
    background: var(--surface-bg, rgba(255, 255, 255, 0.95));
    border: 1px solid var(--border);
    border-radius: 4px;
    box-shadow: 0 16px 24px var(--shadow, rgba(0, 0, 0, 0.15));
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
