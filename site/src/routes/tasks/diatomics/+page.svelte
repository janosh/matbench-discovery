<script lang="ts">
  import { afterNavigate } from '$app/navigation'
  import { page } from '$app/state'
  import { MetricsTable, type ModelData } from '$lib'
  import { ALL_METRICS, DIATOMICS_METRICS, METADATA_COLS } from '$lib/labels'
  import { DiatomicCurve } from '$lib/plot'
  import { sort_from_query, sync_url_params, valid_query_param } from '$lib/url-state'
  import type { SortDir } from '$lib/types'
  import DiatomicsNote from './diatomics-note.md'
  import Select from 'svelte-multiselect'
  import { element_data } from 'matterviz/element'
  import type { ChemicalElement, ElementCategory } from 'matterviz/element'
  import { pick_contrast_color, PLOT_COLORS } from 'matterviz'
  import { ELEM_SYMBOLS } from 'matterviz/labels'
  import { untrack } from 'svelte'
  import { SvelteSet } from 'svelte/reactivity'

  let { data } = $props()
  let diatomic_models = $derived(data?.diatomic_models ?? [])
  let diatomic_curves = $derived(data?.diatomic_curves ?? {})
  let reference_names = $derived(data?.reference_names ?? [])
  let errors = $derived(data?.errors ?? {})
  let error_entries = $derived(Object.entries(errors))
  let error_title = $derived(
    error_entries.map(([model_name, error]) => `${model_name}: ${error}`).join(`\n`),
  )

  const homo_nuc_key = `homo-nuclear`
  const elements = element_data as ChemicalElement[]
  const element_by_symbol = new Map(elements.map((element) => [element.symbol, element]))
  const nonmetal_categories = new Set<ElementCategory>([
    `diatomic nonmetal`,
    `polyatomic nonmetal`,
  ])
  type ElementGroup = {
    key: string
    label: string
    title: string
    includes: (element: ChemicalElement) => boolean
  }
  const category_group = (
    key: string,
    label: string,
    category: ElementCategory,
  ): ElementGroup => ({
    key,
    label,
    title: label,
    includes: (element) => element.category === category,
  })
  const element_groups: ElementGroup[] = [
    { key: `all`, label: `All`, title: `Show all elements`, includes: () => true },
    category_group(`alkali`, `Alkali metals`, `alkali metal`),
    category_group(`alkaline_earth`, `Alkaline earth metals`, `alkaline earth metal`),
    category_group(`transition`, `Transition metals`, `transition metal`),
    category_group(`post_transition`, `Post-transition metals`, `post-transition metal`),
    category_group(`metalloid`, `Metalloids`, `metalloid`),
    {
      key: `nonmetal`,
      label: `Nonmetals`,
      title: `Diatomic and polyatomic nonmetals`,
      includes: (element) => nonmetal_categories.has(element.category),
    },
    {
      key: `halogen`,
      label: `Halogens`,
      title: `Group 17 halogen elements`,
      includes: (element) => element.column === 17 && element.row <= 7,
    },
    category_group(`noble_gas`, `Noble gases`, `noble gas`),
    category_group(`lanthanide`, `Lanthanides`, `lanthanide`),
    category_group(`actinide`, `Actinides`, `actinide`),
  ]
  const element_group_keys = new Set(element_groups.map(({ key }) => key))
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
  const default_sort: { column: string; dir: SortDir } = {
    column: DIATOMICS_METRICS.energy_jump.key,
    dir: `asc`,
  }
  let sort = $state({ ...default_sort })

  // Generate list of homo-nuclear diatomic formulas for elements 1-119
  const homo_diatomic_formulas = ELEM_SYMBOLS.map((symbol) => `${symbol}-${symbol}`)

  // Create a Map to store model colors consistently
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

  let plot_height = $state<number | null>(300)
  let clamped_plot_height = $derived(Math.min(500, Math.max(100, plot_height ?? 300)))
  let selected_element_group = $state(`all`)
  let selected_group = $derived(
    element_groups.find((group) => group.key === selected_element_group) ??
      element_groups[0],
  )

  let url_ready = $state(false)
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

  function default_selected_names(): string[] {
    return [
      // DFT references on by default, plus the newest models with curve data
      ...reference_names,
      ...selectable_names
        .filter((model_name) => !reference_names.includes(model_name))
        .slice(0, default_n_models),
    ]
  }

  const options_from_names = (model_names: string[]) =>
    selectable_options.filter((option) => model_names.includes(String(option.value)))

  let selected_model_options = $state(
    untrack(() => options_from_names(default_selected_names())),
  )
  let selected_model_names = $derived(
    selected_model_options.map((option) => String(option.value)),
  )
  const visible_diatomics = new SvelteSet<string>()
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

  function names_from_url(params: URLSearchParams): string[] {
    const model_param = params.get(`models`)
    if (model_param === null) return default_selected_names()
    if (!model_param) return []

    return model_param
      .split(`,`)
      .map((key_or_ref) => model_name_by_key.get(key_or_ref) ?? key_or_ref)
      .filter((model_name) => selectable_names.includes(model_name))
  }

  const selected_param_value = (model_names: string[]): string =>
    selectable_names
      .filter((model_name) => model_names.includes(model_name))
      .map((model_name) => model_key_by_name.get(model_name) ?? model_name)
      .join(`,`)

  function element_group_from_url(params: URLSearchParams): string {
    return valid_query_param(params, `elements`, `all`, element_group_keys)
  }

  afterNavigate(() => {
    selected_model_options = options_from_names(names_from_url(page.url.searchParams))
    selected_element_group = element_group_from_url(page.url.searchParams)
    sort = sort_from_query(page.url.searchParams, default_sort)
    url_ready = true
  })

  $effect(() => {
    if (!url_ready) return

    const selected = selected_param_value(selected_model_names)
    const defaults = selected_param_value(default_selected_names())
    sync_url_params(
      [
        [`models`, selected, defaults],
        [`elements`, selected_element_group, `all`],
        [`sort`, sort.column, default_sort.column],
        [`dir`, sort.dir, default_sort.dir],
      ],
      page.state,
    )
  })

  function curves_for_formula(formula: string) {
    return selected_model_names.flatMap((model) => {
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
  }

  function observe_plot(node: HTMLElement, formula: string) {
    const hide_plot = () => visible_diatomics.delete(formula)
    const Observer = globalThis.IntersectionObserver
    if (!Observer) {
      visible_diatomics.add(formula)
      return { destroy: hide_plot }
    }
    const observer = new Observer(
      ([entry]) => {
        if (!entry) return
        if (entry.isIntersecting) visible_diatomics.add(formula)
        else visible_diatomics.delete(formula)
      },
      { rootMargin: `300px 0px` },
    )
    observer.observe(node)
    return {
      destroy: () => {
        hide_plot()
        observer.disconnect()
      },
    }
  }
</script>

<h1>Diatomics</h1>

<DiatomicsNote />

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
    <input
      class="plot-height-number"
      type="number"
      min="100"
      max="500"
      step="10"
      bind:value={plot_height}
      aria-label="Plot height in pixels"
    />px
  </label>

  <div class="element-groups" aria-label="Element group filter">
    {#each element_groups as { label, title, key } (key)}
      {@const selected = selected_element_group === key}
      <button
        class:selected
        aria-pressed={selected}
        onclick={() => (selected_element_group = key)}
        {title}
      >
        {label}
      </button>
    {/each}
  </div>

  <Select
    class="model-select"
    options={selectable_options}
    bind:selected={selected_model_options}
    placeholder="Select models to plot"
    style="width: 100%; max-width: 100%; min-width: 19rem"
  />
</div>

<div class="grid" style="--plot-height: {clamped_plot_height}px">
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
    align-items: center;
    gap: 1ex;
  }
  .plot-height-number {
    width: 5em;
  }
  .element-groups {
    display: flex;
    flex-wrap: wrap;
    gap: 0.5em;
    justify-content: center;
  }
  .element-groups button {
    --model-color: var(--text-color, currentColor);
  }
  button {
    padding: 0.3em 0.6em;
    border: none;
    border-radius: 4px;
    background: color-mix(in srgb, var(--model-color, currentColor) 12%, transparent);
    color: var(--text-color, currentColor);
    cursor: pointer;
    font-weight: 500;
  }
  button.selected {
    background: color-mix(in srgb, var(--model-color, currentColor) 35%, transparent);
    font-weight: 650;
  }
  .error-summary {
    margin: 1em auto;
    max-width: 80ch;
    padding: 0.75em 1em;
    border: 1px solid var(--danger, #b91c1c);
    border-radius: 4px;
  }
  button:hover:not(.error, :disabled) {
    background: color-mix(in srgb, var(--model-color, currentColor) 22%, transparent);
  }
</style>
