<script lang="ts">
  import { dev } from '$app/environment'
  import { page } from '$app/stores'
  import per_elem_each_errors from '$figs/per-element-each-errors.json'
  import { MODEL_METADATA, PtableInset } from '$lib'
  import TRAINING_SETS from '$root/data/training-sets.yml'
  import pkg from '$site/package.json'
  import Icon from '@iconify/svelte'
  import type { ChemicalElement } from 'elementari'
  import {
    ColorBar,
    ColorScaleSelect,
    PeriodicTable,
    pretty_num,
    TableInset,
  } from 'elementari'
  import { CopyButton, Tooltip } from 'svelte-zoo'

  $: model_key = $page.params.slug
  $: model = MODEL_METADATA.find(
    (model) => model.model_name.toLowerCase().replaceAll(` `, `-`) == model_key,
  )
  let color_scale: string[] = [`Viridis`]
  let active_element: ChemicalElement | null = null

  // TODO make this dynamic (static n_days_ago from time of last site build is misleading)
  function n_days_ago(dateString: string): string {
    return (
      (new Date().getTime() - new Date(dateString).getTime()) /
      (1000 * 60 * 60 * 24)
    ).toLocaleString(`en-US`, { maximumFractionDigits: 0 })
  }

  export const snapshot = {
    capture: () => ({ color_scale }),
    restore: (values) => ({ color_scale } = values),
  }
</script>

{#if model}
  <div class="model-detail">
    <h1 style="font-size: 2.5em;">{model.model_name}</h1>

    <section class="meta-info">
      <div>
        Version: {#if model.repo}
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
      </div>
      <div>
        <Icon icon="ion:ios-calendar" inline />
        <Tooltip text="{n_days_ago(model.date_added)} days ago">
          <span>Added: {model.date_added}</span>
        </Tooltip>
      </div>
      <div>
        <Icon icon="ri:calendar-check-line" inline />
        <Tooltip text="{n_days_ago(model.date_published)} days ago">
          <span>Published: {model.date_published}</span>
        </Tooltip>
      </div>
      <div>
        <Icon icon="eos-icons:neural-network" inline />
        <Tooltip text={model.model_params.toLocaleString()}>
          <span>{pretty_num(model.model_params, `.3~s`)} parameters</span>
        </Tooltip>
      </div>
      {#if model.n_estimators > 1}
        <div>
          <Icon icon="material-symbols:forest" inline />
          <span>Ensemble {model.n_estimators} models</span>
        </div>
      {/if}
      {#if model.missing_preds}
        <div>
          <Icon icon="fluent:missing-metadata-24-regular" inline />
          <span>
            Missing preds: {pretty_num(model.missing_preds, `,.0f`)}
            <small> ({model.missing_percent})</small>
          </span>
        </div>
      {/if}
      {#if model.pypi}
        <code>
          pip install {model.pypi.split(`/`).pop()}
          <!-- TODO add custom CopyButton labels to remove text -->
          <CopyButton />
        </code>
      {/if}
    </section>

    <section class="links">
      <a href={model.repo} target="_blank" rel="noopener noreferrer">
        <Icon icon="octicon:mark-github" inline /> Repo
      </a>
      <a href={model.paper} target="_blank" rel="noopener noreferrer">
        <Icon icon="ion:ios-paper" inline /> Paper
      </a>
      {#if model.url}
        <a href={model.url} target="_blank" rel="noopener noreferrer">
          <Icon icon="ion:ios-globe" inline /> Docs
        </a>
      {/if}
      <a href={model.doi} target="_blank" rel="noopener noreferrer">
        <Icon icon="academicons:doi" inline /> DOI
      </a>
      <a
        href={`${pkg.repository}/blob/-/models/${model.dirname?.split(`/`).pop()}`}
        target="_blank"
        rel="noopener noreferrer"
      >
        <Icon icon="octicon:file-directory" inline /> Files
      </a>
      {#if model.pypi}
        <a href={model.pypi} target="_blank" rel="noopener noreferrer">
          <Icon icon="simple-icons:pypi" inline /> PyPI
        </a>
      {/if}
    </section>

    {#each [[`e-form`, `Formation Energies`], [`each`, `Convex Hull Distance`]] as [which_energy, title]}
      {#await import(`$figs/energy-parity/${which_energy}-parity-${model.model_name.toLowerCase().replaceAll(` `, `-`)}.svelte`) then ParityPlot}
        <!-- negative margin-bottom corrects for display: none plot title -->
        <h3 style="margin-bottom: -2em;">
          DFT vs ML {title}
        </h3>
        <svelte:component this={ParityPlot.default} height="500" />
      {:catch error}
        {#if dev}
          <p>Failed to load plot:</p>
          <pre>{error}</pre>
        {/if}
      {/await}
    {/each}
    {#if $page.params.slug == `mace`}
      {#await import(`./mace.md`) then Mace}
        <Mace.default />
      {/await}
    {/if}

    <ColorScaleSelect bind:selected={color_scale} />

    <h3>Convex hull distance prediction errors projected onto elements</h3>
    <PeriodicTable
      heatmap_values={per_elem_each_errors[model.model_name]}
      color_scale={color_scale[0]}
      bind:active_element
      tile_props={{ precision: `0.2` }}
      show_photo={false}
    >
      <TableInset slot="inset" style="align-content: center;">
        <PtableInset
          element={active_element}
          elem_counts={per_elem_each_errors[model.model_name]}
          show_percent={false}
          unit="<small style='font-weight: lighter;'>eV / atom</small>"
        />
        <ColorBar
          text_side="top"
          color_scale={color_scale[0]}
          tick_labels={5}
          style="width: 85%; margin: 0 2em;"
        />
      </TableInset>
    </PeriodicTable>

    <section class="authors">
      <h2>Authors</h2>
      <ol>
        {#each model.authors as author}
          <li>
            <span>{author.name}</span>
            {#if author.affiliation}<span class="affiliation">({author.affiliation})</span
              >{/if}
            {#if author.email}<a href="mailto:{author.email}">
                <Icon icon="mdi:email" inline />
              </a>{/if}
            {#if author.github}<a
                href={author.github}
                target="_blank"
                rel="noopener noreferrer"
              >
                <Icon icon="simple-icons:github" inline />
              </a>{/if}
            {#if author.orcid}
              <a href={author.orcid} target="_blank" rel="noopener noreferrer">
                <Icon icon="simple-icons:orcid" inline />
              </a>{/if}
          </li>
        {/each}
      </ol>
    </section>

    {#if model.trained_by}
      <section class="trained-by">
        <h2>Trained By</h2>
        <ol>
          {#each model.trained_by as trainer}
            <li>
              <span>{trainer.name}</span>
              {#if trainer.affiliation}<span class="affiliation"
                  >({trainer.affiliation})</span
                >{/if}
              {#if trainer.orcid}<a
                  href={trainer.orcid}
                  target="_blank"
                  rel="noopener noreferrer"
                >
                  <Icon icon="simple-icons:orcid" inline />
                </a>{/if}
              {#if trainer.github}<a
                  href={trainer.github}
                  target="_blank"
                  rel="noopener noreferrer"
                >
                  <Icon icon="simple-icons:github" inline />
                </a>{/if}
            </li>
          {/each}
        </ol>
      </section>
    {/if}

    <section class="model-info">
      <h2>Model Info</h2>
      <ul>
        {#each [[`Model Version`, model.model_version], [`Model Type`, model.model_type], [`Targets`, model.targets], [`Openness`, model.openness], [`Train Task`, model.train_task], [`Test Task`, model.test_task], [`Trained for Benchmark`, model.trained_for_benchmark ? `Yes` : `No`]] as [key, value]}
          <li>
            {key}
            <strong>{value}</strong>
          </li>
        {/each}
      </ul>
    </section>

    {#if model.training_set}
      <h2>Training Set</h2>
      <section class="training-set">
        {#each Array.isArray(model.training_set) ? model.training_set : [model.training_set] as train_set}
          {@const train_set_info =
            typeof train_set == `string` ? TRAINING_SETS[train_set] : train_set}
          {@const { n_structures, url, title, n_materials } = train_set_info}
          {@const pretty_n_mat =
            typeof n_materials == `number` ? pretty_num(n_materials) : n_materials}
          <p>
            <a href={url} target="_blank" rel="noopener noreferrer">{title}</a>:
            <Tooltip text={n_structures.toLocaleString()}>
              <strong slot="trigger">{pretty_num(n_structures)}</strong>
            </Tooltip>
            structures
            {#if n_materials}
              <Tooltip text={n_materials.toLocaleString()}>
                from <strong slot="trigger">{pretty_n_mat}</strong> materials
              </Tooltip>
            {/if}
          </p>
        {/each}
      </section>
    {/if}

    {#if model.notes}
      <section class="notes">
        {#each Object.entries(model.notes) as [key, note]}
          <h2>{key}</h2>
          {#if typeof note === `string`}
            <p>{@html note}</p>
          {:else if Array.isArray(note)}
            <ol>
              {#each note as val}
                <li>{@html val}</li>
              {/each}
            </ol>
          {:else}
            <ul>
              {#each Object.entries(note) as [key, val]}
                <li><strong>{key}:</strong> {val}</li>
              {/each}
            </ul>
          {/if}
        {/each}
      </section>
    {/if}

    {#if model.hyperparams}
      <section class="hyperparams">
        <h2>Hyperparameters</h2>
        <ul>
          {#each Object.entries(model.hyperparams) as [key, value]}
            <li><strong>{key}:</strong> <code>{JSON.stringify(value)}</code></li>
          {/each}
        </ul>
      </section>
    {/if}

    {#if model.requirements}
      <section class="deps">
        <h2>Dependencies</h2>
        <ul>
          {#each Object.entries(model.requirements) as [pkg, version]}
            {@const href = version.startsWith(`http`)
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
{:else}
  <p>Model not found.</p>
{/if}

<style>
  section {
    margin-bottom: 1em;
  }

  h2 {
    margin: 1em auto 0;
    padding: 0;
  }

  h3 {
    margin: 1em 0;
  }

  section:is(.deps, .model-info) ul {
    display: flex;
    flex-wrap: wrap;
    gap: 1em;
    padding: 0;
  }
  section:is(.deps, .model-info) ul li {
    background-color: rgba(255, 255, 255, 0.1);
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

  .meta-info,
  .links {
    display: flex;
    flex-wrap: wrap;
    gap: 3ex;
    place-content: center;
    margin: 2em auto;
  }

  .links a {
    display: inline-flex;
    align-items: center;
    gap: 5px;
    padding: 5px 10px;
    background-color: rgba(255, 255, 255, 0.1);
    border-radius: 5px;
    text-decoration: none;
    color: lightgray;
  }

  li {
    margin: 1ex 0;
  }

  .affiliation {
    font-style: italic;
    color: gray;
  }

  ul li {
    overflow: hidden;
    white-space: nowrap;
    text-overflow: ellipsis;
  }

  /* hide plotly titles */
  div.model-detail :not(section.notes) :global(h3) {
    text-align: center;
    margin: 2em auto 0;
  }
  div.model-detail :global(svg g.infolayer g.g-gtitle) {
    display: none;
  }
</style>
