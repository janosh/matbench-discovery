<script lang="ts">
  import TestSet from '#lib/benchmark/TestSet.svelte'
  import TaskNavigation from '#lib/benchmark/TaskNavigation.svelte'
  import MetricsTable from '#lib/table/MetricsTable.svelte'
  import ModelSelect from '#lib/ModelSelect.svelte'
  import {
    ACTIVE_MODELS,
    has_diatomics_curves,
    task_page_filters,
  } from '#lib/models.svelte.js'
  import { ButtonGroup } from 'svelte-widgets'
  import { CDS_CONFIG, DEFAULT_CDS_CONFIG } from '#lib/combined-scores.svelte.js'
  import {
    DIATOMICS_METRICS,
    METADATA_COLS,
    scatter_axis_label,
    task_page_visible_cols,
  } from '#lib/labels.js'
  import DiatomicCurve from '#lib/plot/DiatomicCurve.svelte'
  import DynamicScatter from '#lib/plot/DynamicScatter.svelte'
  import ScoreWeights from '#lib/ScoreWeights.svelte'
  import { UrlModelSelection } from '#lib/model-selection.svelte.js'
  import DiatomicsNote from './diatomics-note.md'
  import { element_data, element_groups } from 'matterviz/element'
  import { pick_contrast_color, PLOT_COLORS } from 'matterviz/colors'
  import { SvelteSet } from 'svelte/reactivity'
  import type { PageData } from './$types'
  import {
    DFT_REF_COLORS,
    ElementGroupParam,
    make_plot_observer,
  } from './diatomics.svelte.js'

  let { data }: { data: PageData } = $props()
  // resolve keys from the live MODELS so CDS re-weighting re-ranks the defaults
  const find_model = (key: string) => {
    const model = ACTIVE_MODELS.find(({ model_key }) => model_key === key)
    if (!model) throw new Error(`Unknown diatomics model key: ${key}`)
    return model
  }
  let diatomic_models = $derived(data.diatomic_model_keys.map(find_model))
  let { diatomic_curves, reference_names } = $derived(data)
  let model_idx_by_key = $derived(
    Object.fromEntries(diatomic_models.map(({ model_key }, idx) => [model_key, idx])),
  )
  const label_for = (key: string): string =>
    diatomic_models[model_idx_by_key[key] ?? -1]?.model_name ?? key

  const homo_nuc_key = `homo-nuclear`
  const visible_cols = task_page_visible_cols(...Object.values(DIATOMICS_METRICS))
  // default-sort by the combined diatomics score (CDS), best (highest) first
  let plot = $state({
    x: DIATOMICS_METRICS.diatomics_run_time_sec.key,
    y: DIATOMICS_METRICS.diatomics_combined_score.key,
  })

  // cost-vs-fidelity Pareto: sweep wall time (x) vs CDS (y), size = model params,
  // color = training-set size (the two scaling levers)

  const color_for = (key: string): string =>
    DFT_REF_COLORS[key] ??
    PLOT_COLORS[(model_idx_by_key[key] ?? -1) % PLOT_COLORS.length] ??
    `gray`

  const element_group = new ElementGroupParam()

  let selectable_keys = $derived([...reference_names, ...data.diatomic_model_keys])
  let selectable_options = $derived(
    selectable_keys.map((model_key) => {
      const model_color = color_for(model_key)
      const text_color = pick_contrast_color({ background: model_color })
      return {
        label: `${label_for(model_key)}${reference_names.includes(model_key) ? ` (DFT)` : ``}`,
        value: model_key,
        style: {
          selected: `background: ${model_color}; color: ${text_color};`,
          option: ``,
        },
      }
    }),
  )

  // DFT references plus the three highest CDS-ranked models with available curves.
  let default_selected_keys = $derived([
    ...reference_names,
    ...diatomic_models
      .toSorted(
        (left, right) =>
          (right.metrics?.diatomics?.combined_score ?? -Infinity) -
          (left.metrics?.diatomics?.combined_score ?? -Infinity),
      )
      .slice(0, 3)
      .map(({ model_key }) => model_key),
  ])

  const model_selection = new UrlModelSelection(() => ({
    options: selectable_options,
    defaults: default_selected_keys,
  }))
  let selected_model_keys = $derived(model_selection.values)
  const visible_diatomics = new SvelteSet<string>()
  const observe_plot = make_plot_observer(visible_diatomics)
  let diatomics_to_render = $derived(
    // Only render diatomics where at least one model has data
    element_data
      .filter(element_group.group.includes)
      .map(({ symbol }) => `${symbol}-${symbol}`)
      .filter((formula) =>
        selected_model_keys.some(
          (model_key) =>
            diatomic_curves[model_key]?.[homo_nuc_key]?.[formula]?.energies?.length > 0,
        ),
      ),
  )

  const { filters, in_view } = task_page_filters(
    has_diatomics_curves,
    (params) => {
      model_selection.read(params)
      element_group.read(params)
    },
    () => [model_selection.url_entry, element_group.url_entry],
  )

  const curves_for_formula = (formula: string) =>
    selected_model_keys.flatMap((model_key) => {
      const model_curves = diatomic_curves[model_key]
      const curve = model_curves?.[homo_nuc_key]?.[formula]
      if (!curve?.energies.length) return []
      return [
        {
          model_key,
          label: label_for(model_key),
          // DFT references carry per-formula distances; models share one grid
          distances: curve.distances ?? model_curves.distances,
          energies: curve.energies,
          color: color_for(model_key),
          line_width: reference_names.includes(model_key) ? 2.5 : undefined,
        },
      ]
    })
</script>

<h1 id="diatomics">Diatomics</h1>

<p>
  This task tests diatomic potential-energy curves for agreement with DFT, bond geometry,
  and physical consistency across interatomic separations.
</p>
<p>
  Reference data: <a href="#test-set">PBE and r2SCAN diatomic curves</a>. The combined
  diatomics score (CDS, higher is better) combines Accuracy, Geometry, Speed, and
  Physicality. Reference-relative scores use PBE.
</p>

<TaskNavigation />

<h2 id="leaderboard">Leaderboard</h2>
<p>
  <strong>Interpret with caution:</strong> PBE can be unreliable for stretched diatomics,
  and speed is measured on each model's own GPU (shown on hover), not normalized across
  hardware. Models marked * lack curves and are scored on their remaining elements. See
  <a href="#methodology">methodology</a> for reference quality checks and scoring exclusions.
</p>

<section class="full-bleed">
  <MetricsTable
    model_filter={has_diatomics_curves}
    col_filter={(col) => visible_cols[col.key] ?? true}
    default_sort={{ column: DIATOMICS_METRICS.diatomics_combined_score.key, dir: `desc` }}
    {filters}
    repro_task="diatomics"
  />
</section>

<ScoreWeights
  config={CDS_CONFIG}
  default_config={DEFAULT_CDS_CONFIG}
  title_label={DIATOMICS_METRICS.diatomics_combined_score}
/>

<h2 id="model-comparison" style="text-align: center">
  {@html scatter_axis_label(plot.y)} vs {@html scatter_axis_label(plot.x)}
</h2>
<p>
  This defaults to a cost-vs-fidelity Pareto: each model's full diatomic-sweep wall time
  against its CDS, with marker size showing model parameters and color the training-set
  size. Use the axis/color/size selectors to compare any pair of metrics and metadata.
</p>

<DynamicScatter
  models={ACTIVE_MODELS}
  model_filter={in_view}
  bind:x_key={plot.x}
  bind:y_key={plot.y}
  color_key={METADATA_COLS.n_training_materials.key}
  show_pareto_frontier
  style="height: 800px"
/>

<TestSet task="diatomics">
  <p>
    <a
      href="#diatomic-energy-curves"
      onclick={() => {
        model_selection.selected = selectable_options.filter((option) =>
          reference_names.includes(option.value),
        )
        element_group.key = `all`
      }}>Explore DFT reference curves</a
    > using the curve viewer below. Select models there to compare predictions.
  </p>
</TestSet>

<h2 id="diatomic-energy-curves" style="text-align: center">Diatomic Energy Curves</h2>

<div class="controls">
  <ButtonGroup
    label="Element group filter"
    options={element_groups}
    bind:value={element_group.key}
  />

  <ModelSelect options={selectable_options} bind:value={model_selection.selected} />
</div>

<div class="diatomics-grid bleed-1400">
  {#each diatomics_to_render as formula (formula)}
    {@const is_visible = visible_diatomics.has(formula)}
    <div
      class={[`diatomic-plot-shell`, { 'diatomic-plot-placeholder': !is_visible }]}
      {@attach observe_plot(formula)}
    >
      {#if is_visible}
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

<h2 id="methodology">Methodology</h2>
<DiatomicsNote />

<style>
  h1 {
    margin: 0;
  }
  .controls {
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 1em;
    padding: 1em;
  }
</style>
