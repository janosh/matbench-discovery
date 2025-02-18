<script lang="ts">
  import { DiatomicCurve, MODEL_METADATA, type ModelData } from '$lib'
  import { onMount } from 'svelte'

  // Define a sequence of colors to use for models
  const colors = [
    `#2ca02c`, // green
    `#ff7f0e`, // orange
    `#d62728`, // red
    `#1f77b4`, // blue
    `#9467bd`, // purple
    `#8c564b`, // brown
  ]

  // Get models with diatomic data, sorted by date
  const diatomic_models = MODEL_METADATA.filter(
    (model) => model.metrics?.diatomics?.pred_file,
  ).sort((a, b) => new Date(b.date_added).getTime() - new Date(a.date_added).getTime())

  // Create a Map to store model colors consistently
  const model_colors = new Map(
    diatomic_models.map((model, idx) => [model.model_name, colors[idx % colors.length]]),
  )

  // Start with pre-loaded models selected
  let selected_models = $state(new Set<string>())
  let model_data: Record<
    string,
    Record<string, { energies: number[]; forces: number[][] }>
  > = $state({})
  let distances: Record<string, number[]> = $state({})

  let loading = $state(false)
  let error: string | null = $state(null)
  let [plot_width, plot_height] = $state([400, 300])
  let n_models_to_load: number = 5

  // Cache starts empty since we have no initial data
  const model_cache = new Map()

  async function load_initial_models() {
    loading = true
    error = null

    // Get the two most recent models with diatomic data
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
        model_data[model_name] = cached.data
        distances[model_name] = cached.distances
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
      model_cache.set(model_name, {
        data: data[`homo-nuclear`],
        distances: data[`distances`],
      })
      model_data[model_name] = data[`homo-nuclear`]
      distances[model_name] = data[`distances`]
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

  async function toggle_model(model_name: string) {
    if (selected_models.has(model_name)) {
      selected_models.delete(model_name)
    } else {
      const model = diatomic_models.find((m) => m.model_name === model_name)
      if (!model) return

      selected_models.add(model_name)
      await load_model_data(model)
    }
    selected_models = selected_models // trigger reactivity
  }

  let formulas = $derived(Object.keys(Object.values(model_data)[0] ?? {}))
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
    {#each diatomic_models as { model_name }}
      <button
        class:selected={selected_models.has(model_name)}
        onclick={() => toggle_model(model_name)}
        style:--model-color={model_colors.get(model_name)}
      >
        {model_name}
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
    {#each formulas as formula}
      <DiatomicCurve
        {formula}
        curves={[...selected_models].map((model) => ({
          model_key: model,
          distances: distances[model],
          energies: model_data[model][formula].energies,
          color: model_colors.get(model),
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
