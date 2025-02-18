<script lang="ts">
  import {
    HeatmapTable,
    MODEL_METADATA,
    TRAINING_SETS,
    get_metric_rank_order,
    get_pred_file_urls,
  } from '$lib'
  import { pretty_num } from 'elementari'
  import { click_outside } from 'svelte-zoo/actions'
  import { ALL_METRICS, METADATA_COLS } from './metrics'
  import type { DiscoverySet, HeatmapColumn, ModelData } from './types'

  interface Props {
    discovery_set?: DiscoverySet
    model_filter?: (model: ModelData) => boolean
    col_filter?: (col: HeatmapColumn) => boolean
    [key: string]: unknown
  }
  let {
    discovery_set = `unique_prototypes`,
    model_filter = () => true,
    col_filter = () => true,
    ...rest
  }: Props = $props()

  let active_files: { name: string; url: string }[] = $state([])
  let active_model_name = $state(``)
  let pred_file_modal: HTMLDialogElement | null = $state(null)
  let columns: HeatmapColumn[] = $derived(
    [...ALL_METRICS, ...METADATA_COLS]
      .map((col) => ({
        ...col,
        better: col.better ?? get_metric_rank_order(col.label),
        hidden: !col_filter(col),
      }))
      // Ensure Model column comes first
      .sort((col1, _col2) => (col1.label === `Model` ? -1 : 1)),
  )

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
      data_urls[train_set || title] = training_set_info.download_url || ``

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

  // Transform MODEL_METADATA into table data format
  let metrics_data = $derived(
    MODEL_METADATA.filter(
      (model) => model_filter(model) && model.metrics?.discovery?.[discovery_set],
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
          'κ<sub>SRME</sub>': model.metrics?.phonons?.kappa_103?.κ_SRME,
          'Training Set': format_train_set(model.training_set),
          Params: `<span title="${pretty_num(model.model_params, `,`)} trainable model parameters">${pretty_num(model.model_params)}</span>`,
          Targets: targets_str,
          'Date Added': `<span title="${long_date(model.date_added)}">${model.date_added}</span>`,
          Links: {
            paper: {
              url: model.paper || model.doi,
              title: `Read model paper`,
              icon: `📄`,
            },
            repo: { url: model.repo, title: `View source code`, icon: `📦` },
            pr_url: { url: model.pr_url, title: `View pull request`, icon: `🔗` },
            pred_files: { files: get_pred_file_urls(model), name: model.model_name },
          },
        }
      })
      .sort((row1, row2) => (row2.F1 ?? 0) - (row1.F1 ?? 0)),
  )
</script>

<svelte:window
  onkeydown={(event) => {
    if (event.key === `Escape` && pred_file_modal?.open) {
      pred_file_modal.open = false
      event.preventDefault()
    }
  }}
/>

<HeatmapTable data={metrics_data} {columns} {...rest}>
  {#snippet cell({ col, val })}
    {#if col.label === `Links` && val}
      {@const links = val}
      {#each [links.paper, links.repo, links.pr_url] as link}
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
    {:else if [undefined, null].includes(val)}
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
      {#each active_files as file}
        <li>
          <a href={file.url} target="_blank" rel="noopener noreferrer">
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
