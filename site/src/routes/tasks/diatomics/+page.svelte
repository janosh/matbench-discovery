<script lang="ts">
  import { DiatomicCurve, MODEL_METADATA, type ModelData } from '$lib'
  import { elem_symbols } from 'elementari'
  import { onMount } from 'svelte'

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
  // assign each model a different color
  for (const [idx, model] of MODEL_METADATA.entries()) {
    const color = colors[idx % colors.length]
    model.color = color
  }
  const [humu_nuc_key, _hetero_nuc_key] = [`homo-nuclear`, `hetero-nuclear`] as const

  // Generate list of homonuclear diatomic formulas for elements 1-119
  const homo_diatomic_formulas = elem_symbols.map((symbol) => `${symbol}-${symbol}`)

  // Get models with diatomic data, sorted by most recently submitted
  const diatomic_models = MODEL_METADATA.filter(
    (model) => model.metrics?.diatomics?.pred_file,
  ).sort(
    (m1, m2) => new Date(m2.date_added).getTime() - new Date(m1.date_added).getTime(),
  )

  // Start with pre-loaded models selected
  let selected_models = $state(new Set<string>())
  let model_data: Record<
    string,
    {
      distances: number[]
      [humu_nuc_key]: Record<string, { energies: number[]; forces: number[][] }>
    }
  > = $state({})

  let loading = $state(false)
  let error: string | null = $state(null)
  let [plot_width, plot_height] = $state([400, 300])
  let n_models_to_load: number = 5

  // Cache starts empty since we have no initial data
  const model_cache = new Map()

  async function load_initial_models() {
    loading = true
    error = null

    // Get the two most recently submitted models with diatomic data
    const initial_models = diatomic_models.slice(0, n_models_to_load)

    for (const model of initial_models) {
      try {
        await load_model_data(model)
        selected_models.add(model.model_name)
        selected_models = selected_models // trigger reactivity
      } catch (err) {
        console.error(`Failed to load ${model.model_name}: ${err}`)
        // Continue with next model
      }
    }

    loading = false
  }

  onMount(load_initial_models)

  async function load_model_data(model: ModelData) {
    const { model_name } = model

    try {
      // Check if data is already in cache
      if (model_cache.has(model_name)) {
        const cached = model_cache.get(model_name)!
        model_data[model_name] = cached
        return
      }

      loading = true
      const diatomics = model.metrics?.diatomics
      if (!diatomics || typeof diatomics?.pred_file !== `string`) {
        throw new Error(`No diatomic data found for ${model_name}`)
      }

      let response = await fetch(`/${diatomics.pred_file}`, {
        headers: { 'Accept-Encoding': `gzip` },
      })

      // If local file doesn't exist, try to download from Figshare
      if (!response.ok && diatomics?.pred_file_url) {
        response = await fetch(
          `/api/figshare?url=${encodeURIComponent(diatomics.pred_file_url)}`,
        )
      }

      if (!response.ok)
        throw new Error(
          `HTTP error! status: ${response.status} ${response.statusText} ${response.url}`,
        )

      const data = await response.json()

      // Store in cache and current state
      model_cache.set(model_name, data)
      model_data[model_name] = data
    } catch (err) {
      console.error(`Error loading ${model_name} data:`, err)
      error = err instanceof Error ? err.message : `Failed to load diatomic data`
      // Remove from selected if failed to load
      selected_models.delete(model_name)
      selected_models = selected_models // trigger reactivity
    } finally {
      loading = false
    }
  }

  async function toggle_model(model: ModelData) {
    const { model_name } = model
    if (selected_models.has(model_name)) {
      selected_models.delete(model_name)
    } else {
      await load_model_data(model)
      selected_models.add(model_name)
    }
    selected_models = selected_models // trigger reactivity
  }
</script>

<h1>Diatomic Energy Curves</h1>

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
    {#each diatomic_models as model}
      <button
        class:selected={selected_models.has(model.model_name)}
        onclick={() => toggle_model(model)}
        style:--model-color={model.color}
      >
        {model.model_name}
      </button>
    {/each}
  </div>
</div>

{#if loading}
  <p>Loading diatomic data...</p>
{:else if error}
  <p class="error">Error: {error}</p>
{:else}
  <div class="grid" style="--plot-width: {plot_width}px">
    {#each homo_diatomic_formulas as formula}
      <DiatomicCurve
        {formula}
        curves={[...selected_models]
          .filter((model) => {
            const { energies = [] } = model_data[model][humu_nuc_key][formula] ?? {}
            return energies.length > 0
          })
          .map((model) => ({
            model_key: model,
            distances: model_data[model].distances,
            energies: model_data[model][humu_nuc_key][formula].energies,
            color: MODEL_METADATA.find((m) => m.model_name === model)?.color ?? `gray`,
          }))}
        style="height: {plot_height}px"
      />
    {/each}
  </div>
{/if}

<style>
  h1 {
    margin: 0;
  }
  .grid {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(var(--plot-width, 400px), 1fr));
    gap: 15pt 0;
  }
  .error {
    color: red;
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
  button:hover {
    opacity: 0.8;
  }
</style>
