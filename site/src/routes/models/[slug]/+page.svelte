<script lang="ts">
  import { AuthorBrief, DATASETS, ModelRankCard, PtableInset, SelectToggle } from '$lib'
  import {
    discovery_task_tooltips,
    model_type_tooltips,
    openness_tooltips,
    targets_tooltips,
  } from '$lib/metrics'
  import { has_kappa_parity_model } from '$lib/parity/kappa-parity'
  import { EnergyParityPlot, KappaParityPlot } from '$lib/plot'
  import { get_pred_file_urls } from '$lib/models.svelte'
  import pkg from '$site/package.json'
  import type { ChemicalElement, IconName } from 'matterviz'
  import {
    format_num,
    format_relative_time,
    HeatmapTable,
    Icon,
    ColorBar,
  } from 'matterviz'
  import { PeriodicTable, TableInset } from 'matterviz/periodic-table'
  import type { D3InterpolateName } from 'matterviz/colors'
  import { CopyButton } from 'svelte-multiselect'
  import { click_outside, tooltip } from 'svelte-multiselect/attachments'
  import { SvelteSet } from 'svelte/reactivity'
  import { bind_url_params, valid_query_param } from '$lib/url-state.svelte'
  import type { LoadStatus } from '$lib/asset-loader'
  import { per_element_each_errors as per_elem_each_errors } from '$lib/per-element-errors'
  import type { PageData } from './$types'

  type ModelInfoItem = readonly [key: string, value: string, title?: string | null]
  type ExternalLink = { href?: string; icon: IconName; label: string; title: string }

  let { data }: { data: PageData } = $props()

  // per-system MD breakdown: columns shown only when present in the CSV (older
  // pred files lack n_atoms/cost provenance; pressure is NaN for stress-less systems)
  const md_col_defs = [
    { key: `system`, label: `System`, sticky: true },
    { key: `temperature_kelvin`, label: `T (K)`, format: `.0f` },
    { key: `n_atoms`, label: `Atoms`, format: `.0f` },
    { key: `rdf_error`, label: `ΔRDF (%)`, better: `lower`, format: `.1f` },
    { key: `adf_error`, label: `ΔADF (%)`, better: `lower`, format: `.1f` },
    { key: `vdos_error`, label: `ΔvDOS (%)`, better: `lower`, format: `.1f` },
    {
      key: `pressure_mae`,
      label: `P<sub>MAE</sub> (GPa)`,
      better: `lower`,
      format: `.2f`,
    },
    { key: `pressure_error`, label: `ΔP (%)`, better: `lower`, format: `.1f` },
    {
      key: `energy_rmse`,
      label: `ΔE<sub>RMSE</sub> (meV/atom)`,
      better: `lower`,
      format: `.2f`,
    },
    {
      key: `force_rmse`,
      label: `F<sub>RMSE</sub> (meV/Å)`,
      better: `lower`,
      format: `.1f`,
    },
    { key: `run_time_sec`, label: `Time (s)`, better: `lower`, format: `.3~s` },
    { key: `max_gpu_mem_gb`, label: `VRAM (GB)`, better: `lower`, format: `.2f` },
  ] as const
  let md_rows = $derived(data.md_per_system ?? [])
  let md_cols = $derived(
    md_col_defs.filter((col) => md_rows.some((row) => col.key in row)),
  )

  // static: this page has no color-scale picker (see /models/tmi for one that does)
  const color_scale: D3InterpolateName = `interpolateViridis`
  let active_element: ChemicalElement | null = $state(null)
  // energy-parity tab bar: only the active plot is visible; a tab's plot mounts on
  // first activation and stays mounted (hidden) after, so toggling never reloads
  const energy_parity_options = [
    { value: `e-form`, label: `ML vs DFT Formation Energies` },
    { value: `each`, label: `ML vs DFT Convex Hull Distance` },
  ] as const
  type EnergyTab = (typeof energy_parity_options)[number][`value`]
  const default_energy_tab = energy_parity_options[0].value
  const energy_tab_values = new Set<EnergyTab>(
    energy_parity_options.map((option) => option.value),
  )
  let energy_parity_tab = $state<EnergyTab>(default_energy_tab)
  const mounted_energy_tabs = new SvelteSet<EnergyTab>([default_energy_tab])
  // per-tab load status, reported by each EnergyParityPlot (drives button spinners)
  let energy_parity_statuses = $state<Partial<Record<EnergyTab, LoadStatus>>>({})
  $effect(() => {
    mounted_energy_tabs.add(energy_parity_tab)
  })
  bind_url_params(
    (params) => {
      energy_parity_tab = valid_query_param(
        params,
        `energy_tab`,
        default_energy_tab,
        energy_tab_values,
      )
    },
    () => [[`energy_tab`, energy_parity_tab, default_energy_tab]],
  )
  let { model } = $derived(data)
  let added_ago = $derived(format_relative_time(model.date_added))
  let published_ago = $derived(format_relative_time(model.date_published))
  // rendered in order; links whose href isn't an http(s) URL are skipped
  let external_links: ExternalLink[] = $derived([
    {
      href: model.repo,
      icon: `GitHub`,
      label: `Repo`,
      title: `View source code repository`,
    },
    { href: model.paper, icon: `Paper`, label: `Paper`, title: `Read model paper` },
    { href: model.url, icon: `Docs`, label: `Docs`, title: `View model documentation` },
    { href: model.doi, icon: `DOI`, label: `DOI`, title: `Digital Object Identifier` },
    {
      href: model.dirname
        ? `${pkg.repository}/tree/HEAD/models/${model.dirname}`
        : undefined,
      icon: `Directory`,
      label: `Files`,
      title: `Browse model submission files`,
    },
    { href: model.pypi, icon: `PyPI`, label: `PyPI`, title: `Python package on PyPI` },
    { href: model.pr_url, icon: `PullRequest`, label: `PR`, title: `View pull request` },
    {
      href: model.checkpoint_url,
      icon: `Download`,
      label: `Checkpoint`,
      title: `Download model checkpoint`,
    },
  ])
  let model_info_items: ModelInfoItem[] = $derived([
    [`Model Version`, model.model_version],
    [`Model Type`, model.model_type, model_type_tooltips[model.model_type]],
    [`Targets`, model.targets, targets_tooltips[model.targets]],
    [`Openness`, model.openness, openness_tooltips[model.openness]],
    [`Train Task`, model.train_task, discovery_task_tooltips[model.train_task]],
    [`Test Task`, model.test_task, discovery_task_tooltips[model.test_task]],
    [`Trained for Benchmark`, model.trained_for_benchmark ? `Yes` : `No`],
  ])

  let missing_preds = $derived(
    typeof model.metrics?.discovery === `object`
      ? model.metrics.discovery.unique_prototypes?.missing_preds
      : undefined,
  )
</script>

<div class="model-detail">
  <h1 style="font-size: 2.5em; margin: 0">{model.model_name}</h1>

  <section class="meta-info">
    <span>
      <Icon icon="Versions" />
      Version: {#if model.repo?.startsWith(`http`)}
        <a
          href="{model.repo}/releases/tag/{model.model_version}"
          target="_blank"
          rel="noopener noreferrer"
        >
          {model.model_version}
        </a>
      {:else}
        {model.model_version}
      {/if}
    </span>

    <span title={added_ago} {@attach tooltip()}
      ><Icon icon="Calendar" />
      Added: {model.date_added}
    </span>

    <span title={published_ago} {@attach tooltip()}>
      <Icon icon="CalendarCheck" /> Published: {model.date_published}
    </span>

    <span title={model.model_params.toLocaleString()} {@attach tooltip()}>
      <Icon icon="NeuralNetwork" />
      {format_num(model.model_params, `.3~s`)}
      parameters
    </span>

    {#if model.n_estimators > 1}
      <span><Icon icon="Forest" /> Ensemble of {model.n_estimators} models</span>
    {/if}

    {#if missing_preds != undefined}
      <span
        {@attach tooltip({
          content: `Out of ${format_num(DATASETS.WBM.n_structures, `,`)} WBM structures, ${format_num(missing_preds, `,`)} are missing predictions. This refers only to the discovery task of predicting WBM convex hull distances.`,
        })}
      >
        <Icon icon="MissingMetadata" />
        Missing preds: {format_num(missing_preds, `,.0f`)}
        {#if missing_preds != 0}
          <small>
            ({format_num(missing_preds / DATASETS.WBM.n_structures, `.3~%`)})
          </small>
        {/if}
      </span>
    {/if}

    {#if model.pypi}
      {@const pip_cmd = `pip install ${model.pypi.split(`/`).pop()}`}
      <code style="padding: 0 4pt; place-content: center">
        {pip_cmd}
        <CopyButton
          content={pip_cmd}
          labels={{
            ready: { icon: `Copy`, text: `` },
            success: { icon: `Check`, text: `` },
            error: { icon: `Alert`, text: `` },
          }}
        />
      </code>
    {/if}
  </section>

  <section class="links" {@attach tooltip()}>
    {#each external_links as { href, icon, label, title } (label)}
      {#if href?.startsWith(`http`)}
        <a {href} target="_blank" rel="noopener noreferrer" {title}>
          <Icon {icon} />
          {label}
        </a>
      {/if}
    {/each}
    {#if model.metrics}
      {@const pred_files = get_pred_file_urls(model)}
      {#if pred_files.length > 0}
        <details
          class="pred-files"
          {@attach click_outside({ callback: (node) => (node.open = false) })}
        >
          <summary>
            <Icon icon="Graph" /> Predictions
          </summary>
          <div class="dropdown">
            {#each pred_files as { name, url } (url)}
              <a href={url} target="_blank" rel="noopener noreferrer">
                {@html name}
              </a>
            {/each}
          </div>
        </details>
      {/if}
    {/if}
  </section>

  {#if model.model_key}
    <ModelRankCard model_key={model.model_key} />
  {/if}

  <!-- segmented tab bar doubles as the plot title; the active button shows a
  spinner while its plot's data is still loading -->
  <SelectToggle
    class="energy-parity-tabs"
    bind:selected={energy_parity_tab}
    options={energy_parity_options.map((option) => ({
      ...option,
      loading:
        energy_parity_tab === option.value &&
        energy_parity_statuses[option.value] !== `ready` &&
        energy_parity_statuses[option.value] !== `error`,
    }))}
  />
  <!-- only the default tab's plot mounts on page load; the other mounts on first
  activation and then stays mounted-but-hidden so toggling back is instant (asset
  loads are also promise-cached, and keeping the component alive preserves zoom) -->
  {#each energy_parity_options as { value: energy_kind } (energy_kind)}
    {#if mounted_energy_tabs.has(energy_kind)}
      <EnergyParityPlot
        hidden={energy_parity_tab !== energy_kind}
        {model}
        {energy_kind}
        onstatus={(status) => (energy_parity_statuses[energy_kind] = status)}
      />
    {/if}
  {/each}

  {#if has_kappa_parity_model(model.model_key)}
    <KappaParityPlot {model} />
  {/if}

  {#if md_rows.length > 0}
    <section class="md-per-system">
      <h2 style="text-align: center">Molecular dynamics: per-system breakdown</h2>
      <p>
        Errors of this model's NVT rollouts against each DynaMat v1.0 ab-initio reference
        trajectory (see the <a href="/tasks/md">MD task page</a> for metric definitions). Reveals
        which chemistries a model struggles with, which the leaderboard's cross-system means
        hide.
      </p>
      <HeatmapTable
        data={md_rows}
        columns={md_cols}
        initial_sort="system"
        default_num_format=".3~f"
      />
    </section>
  {/if}

  {#if model.model_key && model.model_key in per_elem_each_errors}
    {@const raw_heatmap = per_elem_each_errors[model.model_key]}
    {@const heatmap_values = Object.fromEntries(
      Object.entries(raw_heatmap).filter(
        (entry): entry is [string, number] => entry[1] !== null,
      ),
    )}
    <h2 style="margin: 1em auto; text-align: center" class="toc-exclude">
      Convex hull distance prediction errors projected onto elements
    </h2>
    <PeriodicTable
      {heatmap_values}
      {color_scale}
      bind:active_element
      tile_props={{ float_fmt: `.2f` }}
      show_photo={false}
      missing={{ color: `rgba(255,255,255,0.3)` }}
    >
      {#snippet inset()}
        <TableInset style="align-content: center">
          <div style="height: 2em">
            {#if active_element}
              <PtableInset
                element={active_element}
                elem_counts={heatmap_values}
                show_percent={false}
                unit="<small style='font-weight: lighter;'>eV / atom</small>"
              />
            {/if}
          </div>
          <ColorBar
            title="|E<sub>ML,hull</sub> - E<sub>DFT,hull</sub>| (eV / atom)"
            title_side="top"
            {color_scale}
            range={[0, Math.max(0, ...Object.values(heatmap_values))]}
            style="width: 80%; margin: 0 2em"
          />
        </TableInset>
      {/snippet}
    </PeriodicTable>
  {/if}

  <section class="authors">
    <h2>Model Authors</h2>
    <ol>
      {#each model.authors as author (author.name)}
        <li>
          <AuthorBrief {author} show_affiliation />
        </li>
      {/each}
    </ol>
  </section>

  {#if model.trained_by}
    <section class="trained-by">
      <h2>Trained By</h2>
      <ol>
        {#each model.trained_by as author (author.name)}
          <li>
            <AuthorBrief {author} show_affiliation />
          </li>
        {/each}
      </ol>
    </section>
  {/if}

  <section class="model-info">
    <h2>Model Info</h2>
    <ul>
      {#each model_info_items as [key, value, title = null] (key)}
        <li {title} {@attach tooltip()}>
          {key}
          {#if key === `Targets`}
            <strong>{@html value.replace(/_(.)/g, `<sub>$1</sub>`)}</strong>
          {:else}
            <strong>{value}</strong>
          {/if}
        </li>
      {/each}
    </ul>
  </section>

  {#if model.training_set}
    <h2>Training Set</h2>
    <section class="training-set">
      {#each model.training_set as dataset_key (dataset_key)}
        {@const dataset = DATASETS[dataset_key]}
        {#if dataset}
          {@const { n_structures, name, slug, n_materials } = dataset}
          <p>
            <a href="/data/{slug}">{name}</a>:
            <span title={n_structures.toLocaleString()} {@attach tooltip()}>
              <strong>{format_num(n_structures)}</strong>
            </span>
            structures
            {#if typeof n_materials == `number`}
              from <span title={n_materials.toLocaleString()} {@attach tooltip()}>
                <strong>{format_num(n_materials)}</strong>
              </span> materials
            {/if}
          </p>
        {:else}
          <p>{dataset_key} (unknown dataset)</p>
        {/if}
      {/each}
    </section>
  {/if}

  {#if model.notes?.html}
    <section class="notes">
      {#each Object.entries(model.notes.html) as [key, note] (key)}
        <h2>{key}</h2>
        <p>{@html note}</p>
      {/each}
    </section>
  {/if}

  {#if model.hyperparams}
    <section class="hyperparams">
      <h2>Hyperparameters</h2>
      <ul>
        {#each Object.entries(model.hyperparams) as [key, value] (key)}
          <li><strong>{key}:</strong> <code>{JSON.stringify(value)}</code></li>
        {/each}
      </ul>
    </section>
  {/if}

  {#if model.requirements}
    <section class="deps">
      <h2>Dependencies</h2>
      <ul>
        {#each Object.entries(model.requirements) as [pkg, version] (pkg)}
          {@const href = version?.startsWith(`http`)
            ? version
            : `https://pypi.org/project/${pkg}/${version}`}
          <li>
            {pkg}
            <a {href} target="_blank" rel="noopener noreferrer">{version}</a>
          </li>
        {/each}
      </ul>
    </section>
  {/if}
</div>

<style>
  h2 {
    margin: 2ex auto 0;
  }
  section {
    text-wrap: balance;
  }
  section:is(.deps, .model-info) ul {
    display: flex;
    flex-wrap: wrap;
    gap: 1em;
    padding: 0;
  }
  section:is(.deps, .model-info) ul li {
    background-color: var(--chip-bg);
    padding: 2pt 6pt;
    border-radius: 3pt;
    text-align: center;
    margin: 0;
    font-weight: lighter;
    max-width: 12em;
  }
  section:is(.deps, .model-info) ul li :is(a, strong) {
    display: block;
    font-weight: bold;
  }
  /* segmented control: buttons fused into one bar with rounded outer corners,
  compact height, sitting directly above the plot as its title */
  :global(.energy-parity-tabs.selection-toggle) {
    gap: 0;
    margin: 2em auto 0.5em;
  }
  :global(.energy-parity-tabs.selection-toggle button) {
    padding: 2px 12px;
    border-radius: 0;
    border-width: 0.5px; /* hairline on HiDPI, incl. the active colored border */
  }
  /* fuse adjacent borders; the active button sits on top so its colored border
  wins the shared edge regardless of which side is selected */
  :global(.energy-parity-tabs.selection-toggle button + button) {
    margin-left: -0.5px;
  }
  :global(.energy-parity-tabs.selection-toggle button.active) {
    position: relative;
    z-index: 1;
  }
  :global(.energy-parity-tabs.selection-toggle button:first-child) {
    border-radius: 9999px 0 0 9999px;
  }
  :global(.energy-parity-tabs.selection-toggle button:last-child) {
    border-radius: 0 9999px 9999px 0;
  }
  /* version numbers as light code, less prominent than the package name */
  section.deps ul li a {
    font-weight: normal;
    font-family: var(--font-mono, monospace);
    font-size: 0.9em;
    opacity: 0.85;
  }
  :is(.meta-info, .links) {
    display: flex;
    flex-wrap: wrap;
    gap: 2ex;
    place-content: center;
    margin: 2em auto;
  }
  .links :is(a, summary) {
    display: inline-flex;
    place-items: center;
    gap: 5px;
    padding: 0 5pt;
    background-color: var(--chip-bg);
    border-radius: 5px;
  }
  .links details {
    position: relative;
    cursor: pointer;
  }
  .links .dropdown {
    position: absolute;
    background-color: var(--page-bg);
    border: 1px solid var(--border);
    border-radius: 5px;
    z-index: 3;
    min-width: max-content;
    box-shadow: 0 0 10px var(--shadow);
  }
  .links .dropdown a {
    display: block;
  }
  li {
    margin: 1ex 0;
  }
  ul li {
    overflow: hidden;
    white-space: nowrap;
    text-overflow: ellipsis;
  }
  div.model-detail :not(section.notes) :global(h3) {
    text-align: center;
    margin: 2em auto 0;
  }
</style>
