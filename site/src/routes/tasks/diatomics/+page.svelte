<script lang="ts">
  import { MetricsTable, ModelSelect, type ModelData, MODELS, SelectToggle } from '$lib'
  import {
    CDS_CONFIG,
    DEFAULT_CDS_CONFIG,
    update_models_cds,
  } from '$lib/combined-scores.svelte'
  import {
    DIATOMICS_METRICS,
    METADATA_COLS,
    scatter_axis_label,
    scatter_options_by_key,
    task_page_visible_cols,
  } from '$lib/labels'
  import { DiatomicCurve, DynamicScatter, RadarChart } from '$lib/plot'
  import { UrlModelSelection } from '$lib/model-selection.svelte'
  import {
    apply_weights_param,
    bind_url_params,
    sort_from_query,
    sort_url_entries,
    valid_query_param,
    weights_to_param,
  } from '$lib/url-state.svelte'
  import type { SortDir } from '$lib/types'
  import DiatomicsNote from './diatomics-note.md'
  import type { ChemicalElement } from 'matterviz/element'
  import { pick_contrast_color, PLOT_COLORS } from 'matterviz'
  import { ELEM_SYMBOLS } from 'matterviz/labels'
  import { SvelteSet } from 'svelte/reactivity'
  import type { PageData } from './$types'
  import { element_by_symbol, element_group_keys, element_groups } from './element-groups'
  import { make_plot_observer } from './observe-plot'

  let { data }: { data: PageData } = $props()
  let diatomic_models = $derived(data?.diatomic_models ?? [])
  let diatomic_curves = $derived(data?.diatomic_curves ?? {})
  let reference_names = $derived(data?.reference_names ?? [])
  let errors = $derived(data?.errors ?? {})
  let error_entries = $derived(Object.entries(errors))
  let error_title = $derived(
    error_entries.map(([model_name, error]) => `${model_name}: ${error}`).join(`\n`),
  )

  const homo_nuc_key = `homo-nuclear`
  const visible_cols = task_page_visible_cols(...Object.values(DIATOMICS_METRICS))
  const has_diatomics_metrics = (model: ModelData): boolean =>
    model.metrics?.diatomics != null
  // default-sort by the combined diatomics score (CDS), best (highest) first
  const default_sort: { column: string; dir: SortDir } = {
    column: DIATOMICS_METRICS.diatomics_combined_score.key,
    dir: `desc`,
  }
  let sort = $state({ ...default_sort })

  // cost-vs-fidelity Pareto: sweep wall time (x) vs CDS (y), size = model params,
  // color = training-set size (the two scaling levers)
  const default_scatter_x = DIATOMICS_METRICS.diatomics_run_time_sec.key
  const default_scatter_y = DIATOMICS_METRICS.diatomics_combined_score.key
  let scatter_x = $state(default_scatter_x)
  let scatter_y = $state(default_scatter_y)

  const homo_diatomic_formulas = ELEM_SYMBOLS.map((symbol) => `${symbol}-${symbol}`)

  let model_colors = $derived(
    new Map<string, string>(
      diatomic_models.map((model, idx) => [
        model.model_name,
        PLOT_COLORS[idx % PLOT_COLORS.length],
      ]),
    ),
  )

  // DFT references get fixed, high-contrast colors (not in the model palette) so they
  // read as ground truth
  const ref_colors: Record<string, string> = { PBE: `#000000`, r2SCAN: `#f032e6` }
  const color_for = (key: string): string =>
    ref_colors[key] ?? model_colors.get(key) ?? `gray`

  let selected_element_group = $state(`all`)
  let selected_group = $derived(
    element_groups.find((group) => group.value === selected_element_group) ??
      element_groups[0],
  )

  // default to the 5 most recently added models with curve data (diatomic_models is
  // sorted newest-first); every other model stays toggleable in the buttons below
  const default_n_models = 5
  const model_name_by_key = $derived(
    new Map(diatomic_models.map((model) => [model.model_key, model.model_name])),
  )
  const model_key_by_name = $derived(
    new Map(diatomic_models.map((model) => [model.model_name, model.model_key])),
  )
  let selectable_names = $derived([
    ...reference_names,
    ...diatomic_models
      .map((model: ModelData) => model.model_name)
      .filter((model_name) => model_name in diatomic_curves && !errors[model_name]),
  ])
  let selectable_options = $derived(
    selectable_names.map((model_name) => {
      const model_color = color_for(model_name)
      const text_color = pick_contrast_color({ bg_color: model_color })
      return {
        label: `${model_name}${reference_names.includes(model_name) ? ` (DFT)` : ``}`,
        value: model_name,
        style: {
          selected: `background: ${model_color}; color: ${text_color};`,
          option: ``,
        },
      }
    }),
  )

  // DFT references on by default, plus the newest models with curve data
  let default_selected_names = $derived([
    ...reference_names,
    ...selectable_names
      .filter((model_name) => !reference_names.includes(model_name))
      .slice(0, default_n_models),
  ])

  const model_selection = new UrlModelSelection(() => ({
    options: selectable_options,
    defaults: default_selected_names,
    // accept model keys (canonical) or display names in the URL, drop unknowns
    from_url: (token) => {
      const model_name = model_name_by_key.get(token) ?? token
      return selectable_names.includes(model_name) ? model_name : undefined
    },
    // serialize display names back to model keys (DFT references have no key)
    to_url: (model_name) => model_key_by_name.get(model_name) ?? model_name,
  }))
  let selected_model_names = $derived(model_selection.values)
  const visible_diatomics = new SvelteSet<string>()
  const observe_plot = make_plot_observer(visible_diatomics)
  let diatomics_to_render = $derived(
    // Only render diatomics where at least one model has data
    homo_diatomic_formulas.filter((formula) => {
      const element_symbol = formula.split(`-`, 1)[0] as ChemicalElement[`symbol`]
      const element = element_by_symbol.get(element_symbol)
      return (
        element &&
        selected_group.includes(element) &&
        selected_model_names.some(
          (model) =>
            diatomic_curves[model]?.[homo_nuc_key]?.[formula]?.energies?.length > 0,
        )
      )
    }),
  )

  const read_url_params = (params: URLSearchParams) => {
    model_selection.read(params)
    selected_element_group = valid_query_param(
      params,
      `elements`,
      `all`,
      element_group_keys,
    )
    sort = sort_from_query(params, default_sort)
    scatter_x = valid_query_param(params, `x`, default_scatter_x, scatter_options_by_key)
    scatter_y = valid_query_param(params, `y`, default_scatter_y, scatter_options_by_key)
    apply_weights_param(params.get(`weights`), CDS_CONFIG, DEFAULT_CDS_CONFIG)
  }
  bind_url_params(read_url_params, () => [
    model_selection.url_entry,
    [`elements`, selected_element_group, `all`],
    ...sort_url_entries(sort, default_sort),
    [`x`, scatter_x, default_scatter_x],
    [`y`, scatter_y, default_scatter_y],
    // custom CDS pillar weights (accuracy,geometry,speed,physicality); omitted at defaults
    [`weights`, weights_to_param(CDS_CONFIG, DEFAULT_CDS_CONFIG)],
  ])

  const curves_for_formula = (formula: string) =>
    selected_model_names.flatMap((model) => {
      const model_curves = diatomic_curves[model]
      const curve = model_curves?.[homo_nuc_key]?.[formula]
      if (!curve?.energies.length) return []
      // DFT references carry per-formula distances; models share one grid
      const distances = curve.distances ?? model_curves.distances
      const line_width = reference_names.includes(model) ? 2.5 : undefined
      return [
        {
          model_key: model,
          distances,
          energies: curve.energies,
          color: color_for(model),
          line_width,
        },
      ]
    })
</script>

<h1>Diatomics</h1>

<div class="task-intro">
  <!-- wrapper div: the markdown renders multiple top-level elements which would
  otherwise each become their own flex item -->
  <div><DiatomicsNote /></div>
  <figure class="task-weights">
    <RadarChart
      size={260}
      config={CDS_CONFIG}
      default_config={DEFAULT_CDS_CONFIG}
      title_label={DIATOMICS_METRICS.diatomics_combined_score}
      on_change={(cfg) => update_models_cds(MODELS, cfg as typeof CDS_CONFIG)}
    />
    <figcaption>
      Drag the knob to reweight the four CDS pillars (see &#9432; for definitions); the
      table updates live. Caveats: Accuracy and Geometry score agreement with PBE, which
      is itself unreliable for stretched diatomics, and Speed compares wall times measured
      on heterogeneous hardware.
    </figcaption>
  </figure>
</div>

<MetricsTable
  model_filter={has_diatomics_metrics}
  col_filter={(col) => visible_cols[col.key] ?? true}
  bind:sort
/>

<h2>{@html scatter_axis_label(scatter_y)} vs {@html scatter_axis_label(scatter_x)}</h2>
<p>
  This defaults to a cost-vs-fidelity Pareto: each model's full diatomic-sweep wall time
  against its CDS, with marker size showing model parameters and color the training-set
  size. Use the axis/color/size selectors to compare any pair of metrics and metadata.
</p>

<DynamicScatter
  models={MODELS}
  model_filter={has_diatomics_metrics}
  bind:x_key={scatter_x}
  bind:y_key={scatter_y}
  color_key={METADATA_COLS.n_training_materials.key}
  show_pareto_frontier
  style="height: 800px"
/>

<h2>Diatomic Energy Curves</h2>

{#if error_entries.length > 0}
  <p class="error-summary" role="alert" title={error_title}>
    Failed to load diatomics data for {error_entries.length}
    {error_entries.length === 1 ? `model` : `models`}.
  </p>
{/if}

<div class="controls">
  <SelectToggle
    aria-label="Element group filter"
    options={element_groups}
    bind:selected={selected_element_group}
  />

  <ModelSelect options={selectable_options} bind:selected={model_selection.selected} />
</div>

<div class="diatomics-grid bleed-1400">
  {#each diatomics_to_render as formula (formula)}
    <div
      class="diatomic-plot-shell"
      class:diatomic-plot-placeholder={!visible_diatomics.has(formula)}
      {@attach observe_plot(formula)}
    >
      {#if visible_diatomics.has(formula)}
        <DiatomicCurve
          {formula}
          curves={curves_for_formula(formula)}
          style="height: 300px"
        />
      {:else}
        <h3 class="diatomic-plot-title">{formula}</h3>
      {/if}
    </div>
  {/each}
</div>

<style>
  h1 {
    margin: 0;
  }
  .task-intro {
    margin-bottom: 1em;
  }
  .controls {
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 1em;
    padding: 1em;
  }
  .error-summary {
    margin: 1em auto;
    max-width: 80ch;
    padding: 0.75em 1em;
    border: 1px solid var(--danger, #b91c1c);
    border-radius: 4px;
  }
</style>
