<script lang="ts">
  import { MetricsTable, type ModelData } from '$lib'
  import { ALL_METRICS, DIATOMICS_METRICS, METADATA_COLS } from '$lib/labels'
  import { DiatomicCurve } from '$lib/plot'
  import type { SortDir } from '$lib/types'
  import { ELEM_SYMBOLS } from 'matterviz/labels'
  import { untrack } from 'svelte'
  import { SvelteSet } from 'svelte/reactivity'

  let { data } = $props()
  let diatomic_models = $derived(data?.diatomic_models ?? [])
  let diatomic_curves = $derived(data?.diatomic_curves ?? {})
  let errors = $derived(data?.errors ?? {})
  let error_entries = $derived(Object.entries(errors))
  let error_title = $derived(
    error_entries.map(([model_name, error]) => `${model_name}: ${error}`).join(`\n`),
  )

  // Define a sequence of colors to use for models
  const colors = [
    `#2ca02c`, // Green
    `#ff7f0e`, // Orange
    `#d62728`, // Red
    `#1f77b4`, // Blue
    `#9467bd`, // Purple
    `#8c564b`, // Brown
    `#e377c2`, // Pink
    `#7f7f7f`, // Gray
    `#bcbd22`, // Yellow
    `#17becf`, // Cyan
    `#e6194B`, // Red
    `#3cb44b`, // Green
    `#ffe119`, // Yellow
    `#4363d8`, // Blue
    `#f58231`, // Orange 2
    `#911eb4`, // Purple 2
  ] as const
  type ColorType = (typeof colors)[number]

  const homo_nuc_key = `homo-nuclear`
  const visible_cols: Record<string, boolean> = Object.fromEntries([
    ...Object.values(ALL_METRICS).map((col): [string, boolean] => [col.label, false]),
    ...Object.values(METADATA_COLS).map((col): [string, boolean] => [col.label, true]),
    ...Object.values(DIATOMICS_METRICS).map((col): [string, boolean] => [
      col.label,
      true,
    ]),
  ])
  const has_diatomics_metrics = (model: ModelData): boolean =>
    model.metrics?.diatomics != null && typeof model.metrics.diatomics === `object`
  let sort = $state<{ column: string; dir: SortDir }>({
    column: DIATOMICS_METRICS.energy_jump.key,
    dir: `asc`,
  })

  // Generate list of homo-nuclear diatomic formulas for elements 1-119
  const homo_diatomic_formulas = ELEM_SYMBOLS.map((symbol) => `${symbol}-${symbol}`)

  // Create a Map to store model colors consistently
  let model_colors = $derived(
    new Map<string, ColorType>(
      diatomic_models.map((model: ModelData, idx: number) => [
        model.model_name,
        colors[idx % colors.length],
      ]),
    ),
  )

  let plot_height = $state(300)

  // default to the 5 most recently added models with curve data (diatomic_models is
  // sorted newest-first); every other model stays toggleable in the buttons below
  const default_n_models = 5
  const selected_models = new SvelteSet<string>(
    untrack(() =>
      diatomic_models
        .map((model: ModelData) => model.model_name)
        .filter((model_name: string) => model_name in diatomic_curves)
        .slice(0, default_n_models),
    ),
  )
  let selected_model_names = $derived([...selected_models])
  const visible_diatomics = new SvelteSet<string>()
  let diatomics_to_render = $derived(
    // Only render diatomics where at least one model has data
    homo_diatomic_formulas.filter((formula) =>
      selected_model_names.some(
        (model) =>
          diatomic_curves[model]?.[homo_nuc_key]?.[formula]?.energies?.length > 0,
      ),
    ),
  )

  function toggle_model(model: ModelData) {
    const { model_name } = model
    if (errors[model_name]) return
    if (selected_models.has(model_name)) selected_models.delete(model_name)
    else selected_models.add(model_name)
  }

  function curves_for_formula(formula: string) {
    return selected_model_names.flatMap((model) => {
      const model_curves = diatomic_curves[model]
      const curve = model_curves?.[homo_nuc_key]?.[formula]
      if (!curve?.energies.length) return []
      const { distances } = model_curves
      const color = model_colors.get(model) ?? `gray`
      return [{ model_key: model, distances, energies: curve.energies, color }]
    })
  }

  function observe_plot(node: HTMLElement, formula: string) {
    const Observer = globalThis.IntersectionObserver
    if (!Observer) {
      visible_diatomics.add(formula)
      return
    }
    const observer = new Observer(
      ([entry]) => {
        if (!entry?.isIntersecting) return
        visible_diatomics.add(formula)
        observer.disconnect()
      },
      { rootMargin: `800px` },
    )
    observer.observe(node)
    return {
      destroy: () => observer.disconnect(),
    }
  }
</script>

<h1>Diatomics</h1>

<blockquote>
  Disclaimer: Take the results on this page with a grain of salt. They are still being
  checked for correctness.
</blockquote>

<MetricsTable
  model_filter={has_diatomics_metrics}
  col_filter={(col) => visible_cols[col.label] ?? true}
  bind:sort
/>

<h2>Diatomic Energy Curves</h2>

{#if error_entries.length > 0}
  <p class="error-summary" role="alert" title={error_title}>
    Failed to load diatomics data for {error_entries.length}
    {error_entries.length === 1 ? `model` : `models`}.
  </p>
{/if}

<div class="controls">
  <label class="plot-height-control">
    Plot height:
    <input type="range" min="100" max="500" bind:value={plot_height} />
    {plot_height}px
  </label>

  <div class="model-toggles">
    {#each diatomic_models as model (model.model_name)}
      {@const error = errors[model.model_name]}
      <button
        class:selected={selected_models.has(model.model_name)}
        class:error
        onclick={() => toggle_model(model)}
        disabled={Boolean(error)}
        style:--model-color={model_colors.get(model.model_name) ?? `gray`}
        title={error}
      >
        {model.model_name}
      </button>
    {/each}
  </div>
</div>

<div class="grid" style="--plot-height: {plot_height}px">
  {#each diatomics_to_render as formula (formula)}
    <div class="plot-shell" use:observe_plot={formula}>
      {#if visible_diatomics.has(formula)}
        <DiatomicCurve
          {formula}
          curves={curves_for_formula(formula)}
          style="height: var(--plot-height)"
        />
      {:else}
        <div class="plot-placeholder">
          <h3>{formula}</h3>
        </div>
      {/if}
    </div>
  {/each}
</div>

<style>
  h1 {
    margin: 0;
  }
  .grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(min(100%, 26rem), 1fr));
    gap: 45pt 15pt;
  }
  .plot-shell,
  .plot-placeholder {
    min-height: var(--plot-height, 300px);
  }
  .plot-placeholder {
    display: grid;
    place-items: start center;
    color: var(--text-muted, #777);
  }
  .controls {
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 1em;
    padding: 1em;
  }
  .plot-height-control {
    display: flex;
    gap: 1ex;
  }
  .model-toggles {
    display: flex;
    flex-wrap: wrap;
    gap: 0.5em;
    justify-content: center;
  }
  button {
    padding: 0.3em 0.6em;
    border: 2px solid var(--model-color, gray);
    border-radius: 4px;
    background: transparent;
    color: var(--model-color, currentColor);
    cursor: pointer;
    font-weight: 500;
  }
  button.selected {
    background: var(--model-color, gray);
    color: white;
  }
  button.error {
    opacity: 0.5;
    cursor: not-allowed;
    border-style: dashed;
  }
  .error-summary {
    margin: 1em auto;
    max-width: 80ch;
    padding: 0.75em 1em;
    border: 1px solid var(--danger, #b91c1c);
    border-radius: 4px;
  }
  button:hover:not(.error, :disabled) {
    opacity: 0.8;
  }
</style>
