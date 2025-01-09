<script lang="ts">
  import {
    HeatmapTable,
    MODEL_METADATA,
    TRAINING_SETS,
    get_metric_rank_order,
    model_is_compliant,
  } from '$lib'
  import modeling_tasks from '$pkg/modeling-tasks.yml'
  import { pretty_num } from 'elementari'
  import { METADATA_COLS, METRICS_COLS } from './metrics'
  import type { HeatmapColumn, ModelData } from './types'

  export let discovery_set: `full_test_set` | `most_stable_10k` | `unique_prototypes` =
    `unique_prototypes`
  export let show_non_compliant: boolean = false
  export let show_energy_only: boolean = false
  export let show_metadata: boolean = true
  export let hide_cols: string[] = []
  export let metadata_cols = METADATA_COLS

  let show_pred_files_modal = false
  let active_files: { name: string; path: string }[] = []
  let active_model_name = ``
  let pred_file_modal: HTMLDialogElement | null = null
  let columns: HeatmapColumn[]
  $: columns = [...METRICS_COLS, ...(show_metadata ? metadata_cols : [])].map((col) => ({
    ...col,
    better: col.better ?? get_metric_rank_order(col.label),
    hidden: hide_cols?.includes(col.label) || col.hidden,
  }))

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

  function get_pred_files(model: ModelData) {
    const files: { name: string; path: string }[] = []
    const raw_file_repo_url = `https://github.com/janosh/matbench-discovery/raw/main`

    function find_pred_files(obj: object, parent_key = ``) {
      if (!obj || typeof obj !== `object`) return

      for (const [key, val] of Object.entries(obj)) {
        if (
          [`pred_file`, `pred_file_url`].includes(key) &&
          val &&
          typeof val === `string`
        ) {
          // Get parent key without _pred_file suffix for label lookup
          const pretty_label = modeling_tasks[parent_key]?.label || parent_key
          const path = val.startsWith(`http`) ? val : `${raw_file_repo_url}/${val}`
          files.push({ name: pretty_label, path })
        } else if (typeof val === `object`) {
          find_pred_files(val, key)
        }
      }
    }

    find_pred_files(model.metrics)
    return files
  }

  const long_date = (date: string): string =>
    new Date(date).toLocaleDateString(undefined, {
      weekday: `long`,
      year: `numeric`,
      month: `long`,
      day: `numeric`,
    })

  // Transform MODEL_METADATA into table data format
  $: metrics_data = MODEL_METADATA.filter(
    (model) =>
      (show_energy_only || model.targets != `E`) &&
      (show_non_compliant || model_is_compliant(model)),
  )
    .map((model) => {
      const metrics = model.metrics?.discovery?.[discovery_set]

      const targets = model.targets.replace(/_(.)/g, `<sub>$1</sub>`)
      const targets_str = `<span title="${targets_tooltips[model.targets]}">${targets}</span>`

      // rename metric keys to pretty labels
      return {
        Model: `<a title="Version: ${model.model_version}" href="/models/${model.model_key}">${model.model_name}</a>`,
        F1: metrics?.F1,
        DAF: metrics?.DAF,
        Prec: metrics?.Precision,
        Acc: metrics?.Accuracy,
        TPR: metrics?.TPR,
        TNR: metrics?.TNR,
        MAE: metrics?.MAE,
        RMSE: metrics?.RMSE,
        'R<sup>2</sup>': metrics?.R2,
        'κ<sub>SRME</sub>': model.metrics?.phonons?.κ_SRME,
        'Training Set': format_train_set(model.training_set),
        Params: `<span title="${pretty_num(model.model_params, `,`)} trainable model parameters">${pretty_num(model.model_params)}</span>`,
        Targets: targets_str,
        'Date Added': `<span title="${long_date(model.date_added)}">${model.date_added}</span>`,
        Links: {
          paper: model.paper || model.doi,
          repo: model.repo,
          pr_url: model.pr_url,
          pred_files: { files: get_pred_files(model), name: model.model_name },
        },
      }
    })
    .sort((row1, row2) => (row2.F1 ?? 0) - (row1.F1 ?? 0))
</script>

<svelte:window
  on:keydown={(event) => {
    if (event.key === `Escape` && show_pred_files_modal) {
      show_pred_files_modal = false
      event.preventDefault()
    }
  }}
/>

<HeatmapTable data={metrics_data} {columns} {...$$restProps}>
  <svelte:fragment slot="cell" let:col let:val>
    {#if col.label === `Links` && val}
      {@const links = val}
      {#if links.paper}
        <a href={links.paper} target="_blank" rel="noopener noreferrer" title="Paper"
          >📄</a
        >
      {/if}
      {#if links.repo}
        <a href={links.repo} target="_blank" rel="noopener noreferrer" title="Code">📦</a>
      {/if}
      {#if links.pr_url}
        <a
          href={links.pr_url}
          target="_blank"
          rel="noopener noreferrer"
          title="Pull Request">🔗</a
        >
      {/if}
      {#if links.pred_files}
        <button
          class="pred-files-btn"
          title="Download model prediction files"
          on:click={() => {
            show_pred_files_modal = true
            active_files = links.pred_files.files
            active_model_name = links.pred_files.name
          }}
        >
          📊
        </button>
      {/if}
    {:else}
      {@html val}
    {/if}
  </svelte:fragment>
</HeatmapTable>

<dialog bind:this={pred_file_modal} open={show_pred_files_modal}>
  <div class="modal-content">
    <button
      class="close-btn"
      on:click={() => (show_pred_files_modal = false)}
      title="Close (or click escape)"
    >
      ×
    </button>
    <h3>Download prediction files for {active_model_name}</h3>
    <ol class="pred-files-list">
      {#each active_files as file}
        <li>
          <a
            href={file.path}
            target="_blank"
            rel="noopener noreferrer"
            on:click={() => (show_pred_files_modal = false)}
          >
            {file.name}
          </a>
        </li>
      {/each}
    </ol>
  </div>
</dialog>

<style>
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
