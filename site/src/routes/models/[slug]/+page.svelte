<script lang="ts">
  import { calculate_days_ago, DATASETS, Icon, IconList, PtableInset } from '$lib'
  import {
    discovery_task_tooltips,
    openness_tooltips,
    targets_tooltips,
  } from '$lib/metrics'
  import { get_pred_file_urls } from '$lib/models.svelte'
  import type { ModelData } from '$lib/types'
  import pkg from '$site/package.json'
  import type { ChemicalElement } from 'matterviz'
  import { ColorBar, format_num, PeriodicTable, TableInset } from 'matterviz'
  import type { D3InterpolateName } from 'matterviz/colors'
  import { CopyButton } from 'svelte-multiselect'
  import { click_outside, tooltip } from 'svelte-multiselect/attachments'
  import per_elem_each_errors from '../per-element-each-errors.json'

  interface Props {
    data: { model: ModelData }
  }
  let { data }: Props = $props()

  let color_scale = $state<D3InterpolateName>(`interpolateViridis`)
  let active_element: ChemicalElement | null = $state(null)
  let days_added = calculate_days_ago(data.model.date_added ?? ``)
  let days_published = calculate_days_ago(data.model.date_published ?? ``)

  export const snapshot = {
    capture: () => ({ color_scale }),
    restore: (values) => ({ color_scale } = values),
  }
</script>

{#if data.model}
  {@const model = data.model}
  {@const { missing_preds } = model.metrics?.discovery?.unique_prototypes ??
    {}}
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

      <span title="{days_added} days ago" {@attach tooltip()}><Icon icon="Calendar" />
        Added: {model.date_added}
      </span>

      <span title="{days_published} days ago" {@attach tooltip()}>
        <Icon icon="CalendarCheck" /> Published: {model.date_published}
      </span>

      <span title={model.model_params.toLocaleString()} {@attach tooltip()}>
        <Icon icon="NeuralNetwork" /> {format_num(model.model_params, `.3~s`)}
        parameters
      </span>

      {#if model.n_estimators > 1}
        <span><Icon icon="Forest" /> Ensemble {model.n_estimators} models</span>
      {/if}

      {#if missing_preds != undefined}
        <span
          {@attach tooltip()}
          title="Out of {format_num(DATASETS.WBM.n_structures, `,`)} WBM structures, {format_num(missing_preds, `,`)} are missing predictions. This refers only to the discovery task of predicting WBM convex hull distances."
        >
          <Icon icon="MissingMetadata" />
          Missing preds: {format_num(missing_preds, `,.0f`)}
          {#if missing_preds != 0}
            <small>
              ({format_num(missing_preds / DATASETS.WBM.n_structures, `.3~%`)})</small>
          {/if}
        </span>
      {/if}

      {#if model.pypi}
        <code style="padding: 0 4pt; place-content: center">
          pip install {model.pypi.split(`/`).pop()}
          <CopyButton
            content={`pip install ${model.pypi.split(`/`).pop()}`}
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
      {#if model.repo?.startsWith(`http`)}
        <a
          href={model.repo}
          target="_blank"
          rel="noopener noreferrer"
          title="View source code repository"
        >
          <Icon icon="GitHub" /> Repo
        </a>
      {/if}
      {#if model.paper?.startsWith(`http`)}
        <a
          href={model.paper}
          target="_blank"
          rel="noopener noreferrer"
          title="Read model paper"
        >
          <Icon icon="Paper" /> Paper
        </a>
      {/if}
      {#if model.url?.startsWith(`http`)}
        <a
          href={model.url}
          target="_blank"
          rel="noopener noreferrer"
          title="View model documentation"
        >
          <Icon icon="Docs" /> Docs
        </a>
      {/if}
      {#if model.doi?.startsWith(`http`)}
        <a
          href={model.doi}
          target="_blank"
          rel="noopener noreferrer"
          title="Digital Object Identifier"
        >
          <Icon icon="DOI" /> DOI
        </a>
      {/if}
      <a
        href="{pkg.repository}/blob/-/models/{model.dirname?.split(`/`).pop()}"
        target="_blank"
        rel="noopener noreferrer"
        title="Browse model submission files"
      >
        <Icon icon="Directory" /> Files
      </a>
      {#if model.pypi?.startsWith(`http`)}
        <a
          href={model.pypi}
          target="_blank"
          rel="noopener noreferrer"
          title="Python package on PyPI"
        >
          <Icon icon="PyPI" /> PyPI
        </a>
      {/if}
      <a
        href={model.pr_url}
        target="_blank"
        rel="noopener noreferrer"
        title="View pull request"
      >
        <Icon icon="PullRequest" /> PR
      </a>
      {#if model.checkpoint_url?.startsWith(`http`)}
        <a
          href={model.checkpoint_url}
          target="_blank"
          rel="noopener noreferrer"
          title="Download model checkpoint"
        >
          <Icon icon="Download" /> Checkpoint
        </a>
      {/if}
      {#if model.metrics}
        {@const pred_files = get_pred_file_urls(model)}
        {#if pred_files.length > 0}
          <details
            class="pred-files"
            {@attach click_outside({ callback: (node) => node.open = false })}
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

    <!-- check if Plotly undefined needed for model-page.test.ts since vitest with JSDOM doesn't mock some Browser APIs that Plotly needs -->
    {#if typeof globalThis.Plotly != `undefined`}
      {#each [[`e-form`, `Formation Energies`], [`each`, `Convex Hull Distance`]] as
        [which_energy, title]
        (which_energy)
      }
        {#await import(`$figs/energy-parity/${which_energy}-parity-${model.model_key}.svelte`)
          then ParityPlot
        }
          <!-- negative margin-bottom corrects for display: none plot title -->
          <h3 style="margin: 1em auto -2em; text-align: center">
            DFT vs ML {title}
          </h3>
          <ParityPlot.default height="500" />
        {/await}
      {/each}
    {/if}

    {#if model.model_name in per_elem_each_errors}
      {@const heatmap_values = per_elem_each_errors?.[model.model_name]}
      <h3 style="margin: 1em auto -1em; text-align: center">
        Convex hull distance prediction errors projected onto elements
      </h3>
      <PeriodicTable
        {heatmap_values}
        {color_scale}
        bind:active_element
        tile_props={{ precision: `.2` }}
        show_photo={false}
        missing_color="rgba(255,255,255,0.3)"
      >
        {#snippet inset()}
          {@const style = `height: 2em; visibility: ${active_element ? `visible` : `hidden`};`}
          <TableInset style="align-content: center">
            <PtableInset
              element={active_element}
              elem_counts={heatmap_values}
              show_percent={false}
              unit="<small style='font-weight: lighter;'>eV / atom</small>"
              {style}
            />
            <ColorBar
              title="|E<sub>ML,hull</sub> - E<sub>DFT,hull</sub>| (eV / atom)"
              title_side="top"
              {color_scale}
              range={[0, Math.max(...(Object.values(heatmap_values) as number[]))]}
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
          {@const org_logo = model.org_logos?.find(
          (logo) => logo.name === author.affiliation,
        )}
          <li>
            <span>{author.name}</span>
            {#if author.affiliation}<span class="affiliation">
                &ensp;{author.affiliation}
                {#if org_logo}
                  &nbsp;<IconList icons={[org_logo]} />
                {/if}
              </span>{/if}
            {#if author.email}<a href="mailto:{author.email}" aria-label="Email">
                &nbsp;<Icon icon="Contact" />
              </a>{/if}
            {#if author.github}<a
                href={author.github}
                target="_blank"
                rel="noopener noreferrer"
                aria-label="GitHub"
              >
                &nbsp;<Icon icon="GitHub" />
              </a>{/if}
            {#if author.orcid}
              <a
                href={author.orcid}
                target="_blank"
                rel="noopener noreferrer"
                aria-label="ORCID"
              >
                &nbsp;<Icon icon="Orcid" />
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
                >({trainer.affiliation})</span>{/if}
              {#if trainer.orcid}<a
                  href={trainer.orcid}
                  target="_blank"
                  rel="noopener noreferrer"
                  aria-label="ORCID"
                >
                  <Icon icon="Orcid" />
                </a>{/if}
              {#if trainer.github}<a
                  href={trainer.github}
                  target="_blank"
                  rel="noopener noreferrer"
                  aria-label="GitHub"
                >
                  <Icon icon="GitHub" />
                </a>{/if}
            </li>
          {/each}
        </ol>
      </section>
    {/if}

    <section class="model-info">
      <h2>Model Info</h2>
      <ul>
        {#each [
          [`Model Version`, model.model_version],
          [`Model Type`, model.model_type],
          [`Targets`, model.targets, targets_tooltips[model.targets]],
          [`Openness`, model.openness, openness_tooltips[model.openness]],
          [
            `Train Task`,
            model.train_task,
            discovery_task_tooltips[model.train_task],
          ],
          [
            `Test Task`,
            model.test_task,
            discovery_task_tooltips[model.test_task],
          ],
          [`Trained for Benchmark`, model.trained_for_benchmark ? `Yes` : `No`],
        ] as
          [key, value, title = null]
          (key)
        }
          <li {title} {@attach tooltip()}>
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
    background-color: var(--card-bg);
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
  .meta-info, .links {
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
    background-color: var(--card-bg);
    border-radius: 5px;
  }
  .links details {
    position: relative;
    cursor: pointer;
  }
  .links .dropdown {
    position: absolute;
    margin-top: 5px;
    background-color: var(--page-bg);
    border: 1px solid rgba(255, 255, 255, 0.1);
    border-radius: 5px;
    z-index: 3;
    min-width: max-content;
  }
  .links .dropdown a {
    display: block;
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
