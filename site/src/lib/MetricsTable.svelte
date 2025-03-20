<script lang="ts">
  import {
    HeatmapTable,
    MODEL_METADATA,
    TRAINING_SETS,
    TableControls,
    get_metric_rank_order,
    get_pred_file_urls,
    model_is_compliant,
  } from '$lib'
  import { pretty_num } from 'elementari'
  import { click_outside } from 'svelte-zoo/actions'
  import {
    ALL_METRICS,
    DEFAULT_COMBINED_METRIC_CONFIG,
    METADATA_COLS,
    calculate_combined_score,
  } from './metrics'
  import type {
    CombinedMetricConfig,
    DiscoverySet,
    HeatmapColumn,
    LinkData,
    ModelData,
    TableData,
  } from './types'

  interface Props {
    discovery_set?: DiscoverySet
    model_filter?: (model: ModelData) => boolean
    col_filter?: (col: HeatmapColumn) => boolean
    show_energy_only?: boolean
    show_noncompliant?: boolean
    config?: CombinedMetricConfig
    [key: string]: unknown
  }
  let {
    discovery_set = `unique_prototypes`,
    model_filter = () => true,
    col_filter = () => true,
    show_energy_only = false,
    show_noncompliant = false,
    config = DEFAULT_COMBINED_METRIC_CONFIG,
    ...rest
  }: Props = $props()

  let active_files: { name: string; url: string }[] = $state([])
  let active_model_name = $state(``)
  let pred_file_modal: HTMLDialogElement | null = $state(null)

  // Make metric_config reactive with $state to properly handle updates
  let metric_config = $state({ ...config })

  // Update metric_config when config prop changes
  $effect(() => {
    metric_config = { ...config }
  })

  // Generate tooltip for combined score that shows current weights
  function get_combined_score_tooltip(): string {
    const weights = metric_config.weights.map((w) => w.label).join(`, `)
    return `Combined Performance Score = weighted average of ${weights}`
  }

  // Define CPS column with tooltip
  let combined_score_column: HeatmapColumn = {
    label: `CPS`,
    tooltip: get_combined_score_tooltip(),
    style: `border-right: 1px solid black;`,
    format: `.3f`,
    better: `higher`,
  }

  // Simple approach for column visibility
  let visible_cols = $state<Record<string, boolean>>({})
  let columns = $derived(
    // Use col_filter, combined_score_column, and visible_cols as dependencies
    [combined_score_column, ...ALL_METRICS, ...METADATA_COLS]
      .map((col) => {
        const better = col.better ?? get_metric_rank_order(col.label)

        // append better=higher/lower to tooltip if applicable
        let tooltip = col.tooltip || ``
        if (better === `higher` || better === `lower`) {
          tooltip = tooltip ? `${tooltip} (${better}=better)` : `${better}=better`
        }
        // Use visible_cols if available, otherwise use col_filter
        const hidden =
          col.label in visible_cols ? !visible_cols[col.label] : !col_filter(col)
        return { ...col, better, tooltip, hidden } as HeatmapColumn
      })
      // Ensure Model column comes first
      .sort((col1, _col2) => (col1.label === `Model` ? -1 : 1)),
  )

  // Initialize visible columns
  $effect(() => {
    // Only initialize once when columns are available
    if (columns.length > 0 && Object.keys(visible_cols).length === 0) {
      columns.forEach((col) => {
        visible_cols[col.label] = !col.hidden
      })
    }
  })

  function format_train_set(model_training_sets: string[]) {
    let [total_structs, total_materials] = [0, 0]
    const data_urls: Record<string, string> = {}
    const tooltip: string[] = []

    for (const train_set of model_training_sets) {
      if (!(train_set in TRAINING_SETS)) {
        console.warn(`Training set ${train_set} not found in TRAINING_SETS`)
        continue
      }
      const training_set_info = TRAINING_SETS[train_set]
      const n_structs = training_set_info.n_structures
      const n_materials = training_set_info.n_materials ?? n_structs

      total_structs += n_structs
      total_materials += n_materials

      const title = training_set_info.title || train_set
      data_urls[train_set || title] = training_set_info.url || ``

      if (n_materials !== n_structs) {
        tooltip.push(
          `${title}: ${pretty_num(n_materials, `,`)} materials (${pretty_num(n_structs, `,`)} structures)`,
        )
      } else {
        tooltip.push(`${title}: ${pretty_num(n_materials, `,`)} materials`)
      }
    }

    const data_links = Object.entries(data_urls).map(([key, href]) => {
      if (href) {
        return `<a href="${href}" target="_blank" rel="noopener noreferrer">${key}</a>`
      }
      return key
    })

    const data_str = data_links.join(`+`)
    const new_line = `&#013;` // line break that works in title attribute
    const dataset_tooltip =
      tooltip.length > 1 ? `${new_line}• ${tooltip.join(new_line + `• `)}` : ``

    let title = `${pretty_num(total_materials, `,`)} materials in training set${new_line}${dataset_tooltip}`
    let train_size_str = `<span title="${title}" data-sort-value="${total_materials}">${pretty_num(total_materials)} (${data_str})</span>`

    if (total_materials !== total_structs) {
      title =
        `${pretty_num(total_materials, `,`)} materials in training set ` +
        `(${pretty_num(total_structs, `,`)} structures counting all DFT relaxation ` +
        `frames per material)${dataset_tooltip}`

      train_size_str =
        `<span title="${title}" data-sort-value="${total_materials}">` +
        `${pretty_num(total_materials)} <small>(${pretty_num(total_structs)})</small> ` +
        `(${data_str})</span>`
    }

    return train_size_str
  }

  const targets_tooltips: Record<ModelData[`targets`], string> = {
    E: `Energy`,
    EF_G: `Energy with gradient-based forces`,
    EF_D: `Energy with direct forces`,
    EFS_G: `Energy with gradient-based forces and stress`,
    EFS_D: `Energy with direct forces and stress`,
    EFS_GM: `Energy with gradient-based forces, stress, and magmoms`,
    EFS_DM: `Energy with direct forces, stress, and magmoms`,
  }

  const long_date = (date: string): string =>
    new Date(date).toLocaleDateString(undefined, {
      weekday: `long`,
      year: `numeric`,
      month: `long`,
      day: `numeric`,
    })

  // Helper to safely access nested geo_opt metrics
  function get_geo_opt_property<T>(
    geo_opt: unknown,
    symprec: string,
    property: string,
  ): T | undefined {
    if (!geo_opt || typeof geo_opt !== `object`) return undefined

    const symprec_key = `symprec=${symprec}`
    const metrics = geo_opt[symprec_key as keyof typeof geo_opt]

    if (!metrics || typeof metrics !== `object`) return undefined

    return metrics[property as keyof typeof metrics] as T
  }

  // Check if a model is energy-only (has no force or stress predictions)
  function is_energy_only_model(model: ModelData): boolean {
    return model.targets === `E`
  }

  // Check if a model is noncompliant (doesn't follow submission guidelines)
  function is_noncompliant_model(model: ModelData): boolean {
    return !model_is_compliant(model)
  }

  // Create a combined filter function that respects all filtering conditions
  function create_combined_filter(
    show_energy: boolean,
    show_noncomp: boolean,
  ): (model: ModelData) => boolean {
    return (model: ModelData) => {
      // Apply the user-provided model_filter first
      if (!model_filter(model)) {
        return false
      }

      // Filter energy-only models if not shown
      const is_energy = is_energy_only_model(model)
      if (is_energy && !show_energy) {
        return false
      }

      // Filter noncompliant models if not shown
      const is_noncomp = is_noncompliant_model(model)
      if (is_noncomp && !show_noncomp) {
        return false
      }

      return true
    }
  }

  // Function to calculate metrics_data with combined score - no caching
  function calculate_metrics_data(config: CombinedMetricConfig) {
    // Get the current filter with current state values
    const current_filter = create_combined_filter(show_energy_only, show_noncompliant)

    // Perform the calculation
    return (
      MODEL_METADATA.filter(
        (model) => current_filter(model) && model.metrics?.discovery?.[discovery_set],
      )
        .map((model) => {
          const metrics = model.metrics?.discovery?.[discovery_set]

          // Get RMSD from geo_opt metrics if available, using the first symprec value
          const geo_opt_metrics = model.metrics?.geo_opt
          let rmsd = undefined
          if (geo_opt_metrics && typeof geo_opt_metrics === `object`) {
            // Try to find the first symprec key and get its RMSD
            const symprec_keys = Object.keys(geo_opt_metrics).filter((k) =>
              k.startsWith(`symprec=`),
            )
            if (symprec_keys.length > 0) {
              const symprec_key = symprec_keys[0]
              rmsd = get_geo_opt_property<number>(
                geo_opt_metrics,
                symprec_key.replace(`symprec=`, ``),
                `rmsd`,
              )
            }
          }

          // Get kappa from phonon metrics
          const phonons = model.metrics?.phonons
          const kappa =
            phonons && typeof phonons === `object` && `kappa_103` in phonons
              ? (phonons.kappa_103?.κ_SRME as number | undefined)
              : undefined

          // Calculate combined score
          const combined_score = calculate_combined_score(
            metrics?.F1,
            rmsd,
            kappa,
            config,
          )

          const targets = model.targets.replace(/_(.)/g, `<sub>$1</sub>`)
          const targets_str = `<span title="${targets_tooltips[model.targets]}">${targets}</span>`

          return {
            Model: `<a title="Version: ${model.model_version}" href="/models/${model.model_key}">${model.model_name}</a>`,
            CPS: combined_score,
            F1: metrics?.F1,
            DAF: metrics?.DAF,
            Prec: metrics?.Precision,
            Acc: metrics?.Accuracy,
            TPR: metrics?.TPR,
            TNR: metrics?.TNR,
            MAE: metrics?.MAE,
            RMSE: metrics?.RMSE,
            'R<sup>2</sup>': metrics?.R2,
            'κ<sub>SRME</sub>': kappa,
            RMSD: rmsd,
            'Training Set': format_train_set(model.training_set),
            Params: `<span title="${pretty_num(model.model_params, `,`)} trainable model parameters" data-sort-value="${model.model_params}">${pretty_num(model.model_params)}</span>`,
            Targets: targets_str,
            'Date Added': `<span title="${long_date(model.date_added)}" data-sort-value="${new Date(model.date_added).getTime()}">${model.date_added}</span>`,
            // Add Links as a special property
            Links: {
              paper: {
                url: model.paper || model.doi,
                title: `Read model paper`,
                icon: `📄`,
              },
              repo: { url: model.repo, title: `View source code`, icon: `📦` },
              pr_url: { url: model.pr_url, title: `View pull request`, icon: `🔗` },
              pred_files: { files: get_pred_file_urls(model), name: model.model_name },
            } as LinkData,
          }
        })
        // Sort by combined score (descending)
        .sort((row1, row2) => {
          // Handle NaN values (they should be sorted to the bottom)
          const score1 = row1[`CPS`]
          const score2 = row2[`CPS`]

          if (isNaN(score1) && isNaN(score2)) return 0
          if (isNaN(score1)) return 1
          if (isNaN(score2)) return -1

          return score2 - score1
        })
    )
  }

  // Make metrics_data explicitly reactive with $state
  let metrics_data = $state<TableData>([])

  // Make metrics_data is recalculated whenever filter settings, props, or metric_config changes
  $effect(() => {
    // When metric_config or filters change, recalculate metrics data
    metrics_data = calculate_metrics_data(metric_config)

    // Update the CPS tooltip
    combined_score_column = {
      ...combined_score_column,
      tooltip: get_combined_score_tooltip(),
    }
  })

  // Handle changes to filter options (energy-only and noncompliant models)
  function handle_filter_change(show_energy: boolean, show_noncomp: boolean) {
    // Update the props directly
    show_energy_only = show_energy
    show_noncompliant = show_noncomp

    // Force immediate recalculation
    metrics_data = calculate_metrics_data(metric_config)
  }
</script>

<svelte:window
  onkeydown={(event) => {
    if (event.key === `Escape` && pred_file_modal?.open) {
      pred_file_modal.open = false
      event.preventDefault()
    }
  }}
/>

<HeatmapTable
  data={metrics_data}
  {columns}
  initial_sort_column="CPS"
  initial_sort_direction="desc"
  sort_hint="Click on column headers to sort table rows"
  {...rest}
>
  {#snippet controls()}
    <TableControls
      {show_energy_only}
      {show_noncompliant}
      bind:visible_cols
      on_filter_change={(show_energy, show_noncomp) => {
        handle_filter_change(show_energy, show_noncomp)
      }}
    />
  {/snippet}

  {#snippet cell({ col, val })}
    {#if col.label === `Links` && val && typeof val === `object` && `paper` in val}
      {@const links = val as LinkData}
      {#each [links.paper, links.repo, links.pr_url] as link (link?.title + link?.url)}
        {#if link?.url}
          <a href={link.url} target="_blank" rel="noopener noreferrer" title={link.title}>
            {link.icon}
          </a>
        {/if}
      {/each}
      {#if links.pred_files}
        <button
          class="pred-files-btn"
          title="Download model prediction files"
          onclick={() => {
            if (!pred_file_modal) return
            pred_file_modal.open = true
            active_files = links.pred_files.files
            active_model_name = links.pred_files.name
          }}
        >
          📊
        </button>
      {/if}
    {:else if typeof val === `number` && col.format}
      {pretty_num(val, col.format)}
    {:else if val === undefined || val === null}
      n/a
    {:else}
      {@html val}
    {/if}
  {/snippet}
</HeatmapTable>

<dialog
  bind:this={pred_file_modal}
  use:click_outside={{
    callback: () => {
      if (pred_file_modal?.open) pred_file_modal.open = false
    },
  }}
>
  <div class="modal-content">
    <button
      class="close-btn"
      onclick={() => {
        if (pred_file_modal?.open) pred_file_modal.open = false
      }}
      title="Close (or click escape)"
    >
      ×
    </button>
    <h3>Download prediction files for {active_model_name}</h3>
    <ol class="pred-files-list">
      {#each active_files as { name, url } (name + url)}
        <li>
          <a href={url} target="_blank" rel="noopener noreferrer">
            {name}
          </a>
        </li>
      {/each}
    </ol>
  </div>
</dialog>

<style>
  /* Make the table protrude into page margins */
  :global(.heatmap-table-container) {
    margin-left: calc((100vw - var(--main-width, 60ch)) / -8) !important;
    margin-right: calc((100vw - var(--main-width, 60ch)) / -8) !important;
    width: auto !important;
  }

  dialog {
    visibility: hidden;
    opacity: 0;
    background: var(--light-bg);
    color: var(--text-color);
    border: none;
    border-radius: 5pt;
    padding: 0;
    max-width: min(90vw, 500px);
  }

  dialog[open] {
    visibility: visible;
    opacity: 1;
    z-index: 2;
  }

  .pred-files-btn {
    background: none;
    padding: 0;
  }

  .modal-content {
    padding: 1em;
  }

  .modal-content h3 {
    margin: 0 0 1ex;
  }

  .pred-files-list {
    margin: 0;
    padding: 0 1em;
  }

  .close-btn {
    position: absolute;
    top: 0;
    right: 0;
    background: none;
    cursor: pointer;
    font-size: 24px;
  }
  .close-btn:hover {
    color: var(--link-color);
  }
</style>
