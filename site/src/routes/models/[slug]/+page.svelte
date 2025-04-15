<script lang="ts">
  import { dev } from '$app/environment'
  import per_elem_each_errors from '$figs/per-element-each-errors.json'
  import { calculate_days_ago, DATASETS, get_pred_file_urls, PtableInset } from '$lib'
  import type { ModelData } from '$lib/types'
  import pkg from '$site/package.json'
  import type { ChemicalElement } from 'elementari'
  import {
    ColorBar,
    ColorScaleSelect,
    PeriodicTable,
    pretty_num,
    TableInset,
  } from 'elementari'
  import { CopyButton, Tooltip } from 'svelte-zoo'
  import { click_outside, titles_as_tooltips } from 'svelte-zoo/actions'

  interface Props {
    data: { model: ModelData }
  }

  let { data }: Props = $props()
  let color_scale = $state([`Viridis`])
  let active_element: ChemicalElement | null = $state(null)
  // TODO fix that calculates days ago from site build time, not time of user visiting page
  let days_added = calculate_days_ago(data.model.date_added ?? ``)
  let days_published = calculate_days_ago(data.model.date_published ?? ``)

  export const snapshot = {
    capture: () => ({ color_scale }),
    restore: (values) => ({ color_scale } = values),
  }
</script>

{#if data.model}
  {@const model = data.model}
  {@const { missing_preds, missing_percent } =
    model.metrics?.discovery?.unique_prototypes ?? {}}
  <div class="model-detail">
    <h1 style="font-size: 2.5em;">{model.model_name}</h1>

    <section class="meta-info">
      <div>
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
      </div>

      <div>
        <svg><use href="#icon-calendar"></use></svg>
        Added: <Tooltip text="{days_added} days ago">
          {model.date_added}
        </Tooltip>
      </div>

      <div>
        <svg><use href="#icon-calendar-check"></use></svg>
        Published: <Tooltip text="{days_published} days ago">
          {model.date_published}
        </Tooltip>
      </div>

      <div>
        <svg><use href="#icon-neural-network"></use></svg>
        <Tooltip text={model.model_params.toLocaleString()}>
          {pretty_num(model.model_params, `.3~s`)}
        </Tooltip> parameters
      </div>

      {#if model.n_estimators > 1}
        <div>
          <svg><use href="#icon-forest"></use></svg>
          <span>Ensemble {model.n_estimators} models</span>
        </div>
      {/if}

      {#if missing_preds != undefined}
        <div>
          <svg><use href="#icon-missing-metadata"></use></svg>
          <span>
            Missing preds: {pretty_num(missing_preds, `,.0f`)}
            {#if missing_preds != 0}
              <small> ({missing_percent})</small>
            {/if}
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
      {#if model.repo.startsWith(`http`)}
        <a
          href={model.repo}
          target="_blank"
          rel="noopener noreferrer"
          title="View source code repository"
          use:titles_as_tooltips
        >
          <svg><use href="#icon-github"></use></svg> Repo
        </a>
      {/if}
      {#if model.paper?.startsWith(`http`)}
        <a
          href={model.paper}
          target="_blank"
          rel="noopener noreferrer"
          title="Read model paper"
          use:titles_as_tooltips
        >
          <svg><use href="#icon-paper"></use></svg> Paper
        </a>
      {/if}
      {#if model.url?.startsWith(`http`)}
        <a
          href={model.url}
          target="_blank"
          rel="noopener noreferrer"
          title="View model documentation"
          use:titles_as_tooltips
        >
          <svg><use href="#icon-docs"></use></svg> Docs
        </a>
      {/if}
      {#if model.doi?.startsWith(`http`)}
        <a
          href={model.doi}
          target="_blank"
          rel="noopener noreferrer"
          title="Digital Object Identifier"
          use:titles_as_tooltips
        >
          <svg><use href="#icon-doi"></use></svg> DOI
        </a>
      {/if}
      <a
        href="{pkg.repository}/blob/-/models/{model.dirname?.split(`/`).pop()}"
        target="_blank"
        rel="noopener noreferrer"
        title="Browse model submission files"
        use:titles_as_tooltips
      >
        <svg><use href="#icon-directory"></use></svg> Files
      </a>
      {#if model.pypi?.startsWith(`http`)}
        <a
          href={model.pypi}
          target="_blank"
          rel="noopener noreferrer"
          title="Python package on PyPI"
          use:titles_as_tooltips
        >
          <svg><use href="#icon-pypi"></use></svg> PyPI
        </a>
      {/if}
      <a
        href={model.pr_url}
        target="_blank"
        rel="noopener noreferrer"
        title="View pull request"
        use:titles_as_tooltips
      >
        <svg><use href="#icon-pull-request"></use></svg> PR
      </a>
      {#if model.checkpoint_url?.startsWith(`http`)}
        <a
          href={model.checkpoint_url}
          target="_blank"
          rel="noopener noreferrer"
          title="Download model checkpoint"
          use:titles_as_tooltips
        >
          <svg><use href="#icon-download"></use></svg> Checkpoint
        </a>
      {/if}
      {#if model.metrics}
        {@const pred_files = get_pred_file_urls(model)}
        {#if pred_files.length > 0}
          <details
            class="pred-files"
            use:click_outside={{
              callback: (node) => {
                if (node.open) node.open = false
              },
            }}
          >
            <summary title="Download model prediction files" use:titles_as_tooltips>
              <svg><use href="#icon-graph"></use></svg> Predictions
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

    <!-- check if Plotly undefined needed for model-page.test.ts since vitest with JSDOM doesn't mock some Browser APIs that Plotly needs -->
    {#if typeof globalThis.Plotly != `undefined`}
      {#each [[`e-form`, `Formation Energies`], [`each`, `Convex Hull Distance`]] as [which_energy, title] (which_energy)}
        {#await import(`$figs/energy-parity/${which_energy}-parity-${model.model_key}.svelte`) then ParityPlot}
          <!-- negative margin-bottom corrects for display: none plot title -->
          <h3 style="margin-bottom: -2em;">
            DFT vs ML {title}
          </h3>
          <ParityPlot.default height="500" />
        {:catch error}
          {#if dev}
            <p>Failed to load plot:</p>
            <pre>{error}</pre>
          {/if}
        {/await}
      {/each}
    {/if}
    <ColorScaleSelect bind:selected={color_scale} />

    {#if model.model_name in per_elem_each_errors}
      {@const heatmap_values = per_elem_each_errors?.[model.model_name]}
      <h3>Convex hull distance prediction errors projected onto elements</h3>
      <PeriodicTable
        {heatmap_values}
        color_scale={color_scale[0]}
        bind:active_element
        tile_props={{ precision: `0.2` }}
        show_photo={false}
      >
        {#snippet inset()}
          <TableInset style="align-content: center;">
            <PtableInset
              element={active_element}
              elem_counts={heatmap_values}
              show_percent={false}
              unit="<small style='font-weight: lighter;'>eV / atom</small>"
            />
            <ColorBar
              label_side="top"
              color_scale={color_scale[0]}
              tick_labels={5}
              style="width: 85%; margin: 0 2em;"
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
            <span>{author.name}</span>
            {#if author.affiliation}<span class="affiliation">({author.affiliation})</span
              >{/if}
            {#if author.email}<a href="mailto:{author.email}" aria-label="Email">
                <svg><use href="#icon-mail"></use></svg>
              </a>{/if}
            {#if author.github}<a
                href={author.github}
                target="_blank"
                rel="noopener noreferrer"
                aria-label="GitHub"
              >
                <svg><use href="#icon-github"></use></svg>
              </a>{/if}
            {#if author.orcid}
              <a
                href={author.orcid}
                target="_blank"
                rel="noopener noreferrer"
                aria-label="ORCID"
              >
                <svg><use href="#icon-orcid"></use></svg>
              </a>{/if}
          </li>
        {/each}
      </ol>
    </section>

    {#if model.trained_by}
      <section class="trained-by">
        <h2>Trained By</h2>
        <ol>
          {#each model.trained_by as trainer (trainer.name)}
            <li>
              <span>{trainer.name}</span>
              {#if trainer.affiliation}<span class="affiliation"
                  >({trainer.affiliation})</span
                >{/if}
              {#if trainer.orcid}<a
                  href={trainer.orcid}
                  target="_blank"
                  rel="noopener noreferrer"
                  aria-label="ORCID"
                >
                  <svg><use href="#icon-orcid"></use></svg>
                </a>{/if}
              {#if trainer.github}<a
                  href={trainer.github}
                  target="_blank"
                  rel="noopener noreferrer"
                  aria-label="GitHub"
                >
                  <svg><use href="#icon-github"></use></svg>
                </a>{/if}
            </li>
          {/each}
        </ol>
      </section>
    {/if}

    <section class="model-info">
      <h2>Model Info</h2>
      <ul>
        {#each [[`Model Version`, model.model_version], [`Model Type`, model.model_type], [`Targets`, model.targets], [`Openness`, model.openness], [`Train Task`, model.train_task], [`Test Task`, model.test_task], [`Trained for Benchmark`, model.trained_for_benchmark ? `Yes` : `No`]] as [key, value] (key)}
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
        {#each model.training_set as dataset_key (dataset_key)}
          {@const dataset = DATASETS[dataset_key]}
          {@const { n_structures, title, slug, n_materials } = dataset}
          <p>
            <a href="/data/{slug}">{title}</a>:
            <Tooltip text={n_structures.toLocaleString()}>
              <strong>{pretty_num(n_structures)}</strong>
            </Tooltip>
            structures
            {#if typeof n_materials == `number`}
              from <Tooltip text={n_materials.toLocaleString()}>
                <strong>{pretty_num(n_materials)}</strong>
              </Tooltip> materials
            {/if}
          </p>
        {/each}
      </section>
    {/if}

    {#if model.notes?.html}
      <section class="notes">
        {#each Object.entries(model.notes.html) as [key, note] (key)}
          <h2>{key}</h2>
          {#if typeof note === `string`}
            <p>{@html note}</p>
          {:else if Array.isArray(note)}
            <ol>
              {#each note as val (val)}
                <li>{@html val}</li>
              {/each}
            </ol>
          {:else}
            <ul>
              {#each Object.entries(note) as [key, val] (key)}
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
  .links :is(a, summary) {
    display: inline-flex;
    align-items: center;
    gap: 5px;
    padding: 5px 10px;
    background-color: rgba(255, 255, 255, 0.1);
    border-radius: 5px;
    text-decoration: none;
    color: lightgray;
  }
  .links details {
    position: relative;
    cursor: pointer;
  }
  .links .dropdown {
    position: absolute;
    margin-top: 5px;
    background-color: var(--light-bg);
    border: 1px solid rgba(255, 255, 255, 0.1);
    border-radius: 5px;
    z-index: 3;
    min-width: max-content;
  }
  .links .dropdown a {
    display: block;
    background: none;
  }
  .links .dropdown a:hover {
    background-color: rgba(255, 255, 255, 0.1);
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
  div.model-detail :not(section.notes) :global(h3) {
    text-align: center;
    margin: 2em auto 0;
  }
  /* hide plotly titles */
  div.model-detail :global(svg g.infolayer g.g-gtitle) {
    display: none;
  }
</style>
