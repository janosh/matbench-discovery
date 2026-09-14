<script lang="ts">
  import TestSet from '$lib/benchmark/TestSet.svelte'
  import TaskNavigation from '$lib/benchmark/TaskNavigation.svelte'
  import MetricsTable from '$lib/table/MetricsTable.svelte'
  import ModelSelect from '$lib/ModelSelect.svelte'
  import {
    ACTIVE_MODELS,
    has_diatomics_curves,
    make_table_filters,
  } from '$lib/models.svelte'
  import { ButtonGroup } from 'svelte-widgets'
  import { CDS_CONFIG, DEFAULT_CDS_CONFIG } from '$lib/combined-scores.svelte'
  import {
    DIATOMICS_METRICS,
    METADATA_COLS,
    scatter_axis_label,
    task_page_visible_cols,
  } from '$lib/labels'
  import DiatomicCurve from '$lib/plot/DiatomicCurve.svelte'
  import DynamicScatter from '$lib/plot/DynamicScatter.svelte'
  import RadarChart from '$lib/plot/RadarChart.svelte'
  import { UrlModelSelection } from '$lib/model-selection.svelte'
  import { bind_url_params } from '$lib/url-state.svelte'
  import { valid_query_param } from 'svelte-widgets/url-params'
  import DiatomicsNote from './diatomics-note.md'
  import { element_data } from 'matterviz/element'
  import { pick_contrast_color, PLOT_COLORS } from 'matterviz/colors'
  import { SvelteSet } from 'svelte/reactivity'
  import type { PageData } from './$types'
  import { element_group_keys, element_groups } from './element-groups'
  import { make_plot_observer } from './observe-plot'

  let { data }: { data: PageData } = $props()
  let diatomic_models = $derived(data?.diatomic_models ?? [])
  let diatomic_curves = $derived(data?.diatomic_curves ?? {})
  let reference_names = $derived(data?.reference_names ?? [])
  let model_idx_by_key = $derived(
    Object.fromEntries(diatomic_models.map(({ model_key }, idx) => [model_key, idx])),
  )
  const label_for = (key: string): string =>
    diatomic_models[model_idx_by_key[key] ?? -1]?.model_name ?? key
  let errors = $derived(data?.errors ?? {})
  let error_entries = $derived(Object.entries(errors))

  const homo_nuc_key = `homo-nuclear`
  const visible_cols = task_page_visible_cols(...Object.values(DIATOMICS_METRICS))
  const filters = make_table_filters()
  // default-sort by the combined diatomics score (CDS), best (highest) first
  let plot = $state({
    x: DIATOMICS_METRICS.diatomics_run_time_sec.key,
    y: DIATOMICS_METRICS.diatomics_combined_score.key,
  })

  // cost-vs-fidelity Pareto: sweep wall time (x) vs CDS (y), size = model params,
  // color = training-set size (the two scaling levers)

  // DFT references get fixed, high-contrast colors (not in the model palette) so they
  // read as ground truth
  const ref_colors: Record<string, string> = { PBE: `#000000`, r2SCAN: `#f032e6` }
  const color_for = (key: string): string =>
    ref_colors[key] ??
    PLOT_COLORS[(model_idx_by_key[key] ?? -1) % PLOT_COLORS.length] ??
    `gray`

  let selected_element_group = $state(`all`)
  let selected_group = $derived(
    element_groups.find((group) => group.value === selected_element_group) ??
      element_groups[0],
  )

  let available_models = $derived(
    diatomic_models.filter(
      ({ model_key }) => model_key in diatomic_curves && !errors[model_key],
    ),
  )
  let selectable_keys = $derived([
    ...reference_names,
    ...available_models.map(({ model_key }) => model_key),
  ])
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
    ...available_models
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
      .filter(selected_group.includes)
      .map(({ symbol }) => `${symbol}-${symbol}`)
      .filter((formula) =>
        selected_model_keys.some(
          (model_key) =>
            diatomic_curves[model_key]?.[homo_nuc_key]?.[formula]?.energies?.length > 0,
        ),
      ),
  )

  const read_url_params = (params: URLSearchParams) => {
    model_selection.read(params)
    selected_element_group = valid_query_param(
      params,
      `elements`,
      `all`,
      element_group_keys,
    )
    filters.read(params)
  }
  bind_url_params(read_url_params, () => [
    model_selection.url_entry,
    [`elements`, selected_element_group, `all`],
    ...filters.url_entries,
  ])

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
  and speed compares heterogeneous hardware. Models marked * lack curves and are scored on
  their remaining elements. See <a href="#methodology">methodology</a> for reference quality
  checks and scoring exclusions.
</p>

<section class="full-bleed">
  <MetricsTable
    model_filter={has_diatomics_curves}
    col_filter={(col) => visible_cols[col.key] ?? true}
    default_sort={{ column: DIATOMICS_METRICS.diatomics_combined_score.key, dir: `desc` }}
    {filters}
  />
</section>

<details style="margin-block: 1em">
  <summary>Adjust score weights</summary>
  <figure class="task-weights">
    <RadarChart
      size={260}
      config={CDS_CONFIG}
      default_config={DEFAULT_CDS_CONFIG}
      title_label={DIATOMICS_METRICS.diatomics_combined_score}
    />
    <figcaption>
      Drag the knob to reweight the CDS pillars (see &#9432; for definitions); the table
      and plots update live.
    </figcaption>
  </figure>
</details>

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
  model_filter={has_diatomics_curves}
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
        selected_element_group = `all`
      }}>Explore DFT reference curves</a
    > using the curve viewer below. Select models there to compare predictions.
  </p>
</TestSet>

<h2 id="diatomic-energy-curves" style="text-align: center">Diatomic Energy Curves</h2>

{#if error_entries.length > 0}
  <div class="error-summary" role="alert">
    <p>
      Failed to load diatomics data for {error_entries.length}
      {error_entries.length === 1 ? `model` : `models`}.
    </p>
    <details>
      <summary>Asset details</summary>
      <ul>
        {#each error_entries as [key, error] (key)}
          <li><a href="/models/{key}">{label_for(key)}</a>: {error}</li>
        {/each}
      </ul>
    </details>
  </div>
{/if}

<div class="controls">
  <ButtonGroup
    label="Element group filter"
    options={element_groups}
    bind:selected={selected_element_group}
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
<details>
  <summary>Reference provenance, scoring exclusions, and contributions</summary>
  <DiatomicsNote />
</details>

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
  .error-summary {
    overflow-wrap: anywhere;
    margin: 1em auto;
    max-width: 80ch;
    padding: 0.75em 1em;
    border: 1px solid var(--danger, #b91c1c);
    border-radius: 4px;
  }
</style>
