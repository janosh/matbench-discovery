<script lang="ts">
  import { dev } from '$app/environment'
  import TRAINING_SETS from '$data/training-sets.yml'
  import per_elem_each_errors from '$figs/per-element-each-errors.json'
  import { get_pred_file_urls, PtableInset } from '$lib'
  import pkg from '$site/package.json'
  import type { ChemicalElement } from 'elementari'
  import {
    ColorBar,
    ColorScaleSelect,
    PeriodicTable,
    pretty_num,
    TableInset,
  } from 'elementari'
  import 'iconify-icon'
  import { CopyButton, Tooltip } from 'svelte-zoo'
  import { click_outside, titles_as_tooltips } from 'svelte-zoo/actions'

  interface Props {
    data: { model: ModelData }
  }

  let { data }: Props = $props()
  let color_scale = $state([`Viridis`])
  let active_element: ChemicalElement | null = $state(null)

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

{#if data.model}
  {@const model = data.model}
  {@const { missing_preds, missing_percent } =
    model.metrics?.discovery?.unique_prototypes ?? {}}
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
        <iconify-icon icon="ion:ios-calendar" inline></iconify-icon>
        <Tooltip text="{n_days_ago(model.date_added)} days ago">
          <span>Added: {model.date_added}</span>
        </Tooltip>
      </div>
      <div>
        <iconify-icon icon="ri:calendar-check-line" inline></iconify-icon>
        <Tooltip text="{n_days_ago(model.date_published)} days ago">
          <span>Published: {model.date_published}</span>
        </Tooltip>
      </div>
      <div>
        <iconify-icon icon="eos-icons:neural-network" inline></iconify-icon>
        <Tooltip text="{model.model_params.toLocaleString()} trainable parameters">
          <span>{pretty_num(model.model_params, `.3~s`)} parameters</span>
        </Tooltip>
      </div>
      {#if model.n_estimators > 1}
        <div>
          <iconify-icon icon="material-symbols:forest" inline></iconify-icon>
          <span>Ensemble {model.n_estimators} models</span>
        </div>
      {/if}
      {#if missing_preds != undefined}
        <div>
          <iconify-icon icon="fluent:missing-metadata-24-regular" inline></iconify-icon>
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
      <a
        href={model.repo}
        target="_blank"
        rel="noopener noreferrer"
        title="View source code repository"
        use:titles_as_tooltips
      >
        <iconify-icon icon="octicon:mark-github" inline></iconify-icon> Repo
      </a>
      {#if model.paper?.startsWith(`http`)}
        <a
          href={model.paper}
          target="_blank"
          rel="noopener noreferrer"
          title="Read model paper"
          use:titles_as_tooltips
        >
          <iconify-icon icon="ion:ios-paper" inline></iconify-icon> Paper
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
          <iconify-icon icon="ion:ios-globe" inline></iconify-icon> Docs
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
          <iconify-icon icon="academicons:doi" inline></iconify-icon> DOI
        </a>
      {/if}
      <a
        href={`${pkg.repository}/blob/-/models/${model.dirname?.split(`/`).pop()}`}
        target="_blank"
        rel="noopener noreferrer"
        title="Browse model submission files"
        use:titles_as_tooltips
      >
        <iconify-icon icon="octicon:file-directory" inline></iconify-icon> Files
      </a>
      {#if model.pypi}
        <a
          href={model.pypi}
          target="_blank"
          rel="noopener noreferrer"
          title="Python package on PyPI"
          use:titles_as_tooltips
        >
          <iconify-icon icon="simple-icons:pypi" inline></iconify-icon> PyPI
        </a>
      {/if}
      <a
        href={model.pr_url}
        target="_blank"
        rel="noopener noreferrer"
        title="View pull request"
        use:titles_as_tooltips
      >
        <iconify-icon icon="octicon:git-pull-request" inline></iconify-icon> PR
      </a>
      {#if model.checkpoint_url?.startsWith(`http`)}
        <a
          href={model.checkpoint_url}
          target="_blank"
          rel="noopener noreferrer"
          title="Download model checkpoint"
          use:titles_as_tooltips
        >
          <iconify-icon icon="mdi:download-box" inline></iconify-icon> Checkpoint
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
              <iconify-icon icon="octicon:graph" inline></iconify-icon> Predictions
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
                <iconify-icon icon="mdi:email" inline></iconify-icon>
              </a>{/if}
            {#if author.github}<a
                href={author.github}
                target="_blank"
                rel="noopener noreferrer"
                aria-label="GitHub"
              >
                <iconify-icon icon="simple-icons:github" inline></iconify-icon>
              </a>{/if}
            {#if author.orcid}
              <a
                href={author.orcid}
                target="_blank"
                rel="noopener noreferrer"
                aria-label="ORCID"
              >
                <iconify-icon icon="simple-icons:orcid" inline></iconify-icon>
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
                  <iconify-icon icon="simple-icons:orcid" inline></iconify-icon>
                </a>{/if}
              {#if trainer.github}<a
                  href={trainer.github}
                  target="_blank"
                  rel="noopener noreferrer"
                  aria-label="GitHub"
                >
                  <iconify-icon icon="simple-icons:github" inline></iconify-icon>
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
        {#each Array.isArray(model.training_set) ? model.training_set : [model.training_set] as train_set (train_set)}
          {@const train_set_info =
            typeof train_set == `string` ? TRAINING_SETS[train_set] : train_set}
          {@const { n_structures, url, title, n_materials } = train_set_info}
          {@const pretty_n_mat =
            typeof n_materials == `number` ? pretty_num(n_materials) : n_materials}
          <p>
            <a href={url} target="_blank" rel="noopener noreferrer">{title}</a>:
            <Tooltip text={n_structures.toLocaleString()}>
              {#snippet trigger()}
                <strong>{pretty_num(n_structures)}</strong>
              {/snippet}
            </Tooltip>
            structures
            {#if n_materials}
              <Tooltip text={n_materials.toLocaleString()}>
                from {#snippet trigger()}
                  <strong>{pretty_n_mat}</strong>
                {/snippet} materials
              </Tooltip>
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
