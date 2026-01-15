<script lang="ts">
  import { DiatomicCurve, type ModelData } from '$lib'
  import { ELEM_SYMBOLS } from 'matterviz/labels'
  import { SvelteSet } from 'svelte/reactivity'

  let { data } = $props()
  let diatomic_models = $derived(data?.diatomic_models ?? [])
  let diatomic_curves = $derived(data?.diatomic_curves ?? {})
  let errors = $derived(data?.errors ?? {})

  // Define a sequence of colors to use for models
  const colors = [
    `#2ca02c`, // green
    `#ff7f0e`, // orange
    `#d62728`, // red
    `#1f77b4`, // blue
    `#9467bd`, // purple
    `#8c564b`, // brown
    `#e377c2`, // pink
    `#7f7f7f`, // gray
    `#bcbd22`, // yellow
    `#17becf`, // cyan
    `#e6194B`, // red
    `#3cb44b`, // green
    `#ffe119`, // yellow
    `#4363d8`, // blue
    `#f58231`, // orange 2
    `#911eb4`, // purple 2
  ] as const

  type ColorType = (typeof colors)[number]

  const [humu_nuc_key, _hetero_nuc_key] = [`homo-nuclear`, `hetero-nuclear`] as const

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

  let plot_width = $state(400)
  let plot_height = $state(300)

  // Start with pre-loaded models selected, reset when data changes
  const selected_models = new SvelteSet<string>()
  $effect(() => {
    // Clear and repopulate when diatomic_curves changes (e.g., on page load)
    selected_models.clear()
    for (const key of Object.keys(diatomic_curves)) {
      selected_models.add(key)
    }
  })
  let diatomics_to_render = $derived(
    // only render diatomics where at least one model has data
    homo_diatomic_formulas.filter((formula) =>
      [...selected_models].some(
        (model) =>
          diatomic_curves[model]?.[humu_nuc_key]?.[formula]?.energies?.length > 0,
      )
    ),
  )

  function toggle_model(model: ModelData) {
    const { model_name } = model
    if (selected_models.has(model_name)) selected_models.delete(model_name)
    else selected_models.add(model_name)
  }
</script>

<h1>Diatomics</h1>

<blockquote>
  Disclaimer: Take the results on this page with a grain of salt. They are still being
  checked for correctness.
</blockquote>

<h2>Diatomic Energy Curves</h2>

<div class="controls">
  <div class="plot-controls">
    <label>
      Plot width:
      <input type="range" min="200" max="600" bind:value={plot_width} />
      {plot_width}px
    </label>
    <label>
      Plot height:
      <input type="range" min="100" max="500" bind:value={plot_height} />
      {plot_height}px
    </label>
  </div>

  <div class="model-toggles">
    {#each diatomic_models as model (model.model_name)}
      {@const error = errors[model.model_name]}
      <button
        class:selected={selected_models.has(model.model_name)}
        class:error
        onclick={() => toggle_model(model)}
        style:--model-color={model_colors.get(model.model_name) ?? `gray`}
        title={error}
      >
        {model.model_name}
      </button>
    {/each}
  </div>
</div>

<div class="grid" style="--plot-width: {plot_width}px">
  {#each diatomics_to_render as formula (formula)}
    {@const style = `height: ${plot_height}px`}
    <DiatomicCurve
      {formula}
      curves={[...selected_models]
      .filter((model) => {
        const { energies = [] } = diatomic_curves[model]?.[humu_nuc_key]?.[formula] ??
          {}
        return energies.length > 0
      })
      .map((model) => ({
        model_key: model,
        distances: diatomic_curves[model].distances,
        energies: diatomic_curves[model][humu_nuc_key][formula].energies,
        color: model_colors.get(model) ?? `gray`,
      }))}
      {style}
    />
  {/each}
</div>

<style>
  h1 {
    margin: 0;
  }
  .grid {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(var(--plot-width, 400px), 1fr));
    gap: 15pt 0;
  }
  .controls {
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 1em;
    padding: 1em;
  }
  .plot-controls {
    display: flex;
    gap: 2em;
  }
  .plot-controls label {
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
    transition: all 0.2s;
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
  button:hover:not(.error) {
    opacity: 0.8;
  }
</style>
