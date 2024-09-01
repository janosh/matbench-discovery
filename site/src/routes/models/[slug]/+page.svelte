<script lang="ts">
  import { page } from '$app/stores'
  import TRAINING_SETS from '$root/data/training-sets.yml'
  import pkg from '$site/package.json'
  import Icon from '@iconify/svelte'
  import { pretty_num } from 'elementari'
  import { CopyButton, Tooltip } from 'svelte-zoo'
  import Mace from './mace.svelte'

  export let data

  $: model = data.models.find(
    (model) => model.model_name.toLowerCase().replace(/\s+/g, '-') === $page.params.slug,
  )

  function n_days_ago(dateString: string): string {
    return (
      (new Date().getTime() - new Date(dateString).getTime()) /
      (1000 * 60 * 60 * 24)
    ).toLocaleString('en-US', { maximumFractionDigits: 0 })
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
        <Tooltip text={`${n_days_ago(model.date_added)} days ago`}>
          <span>Added: {model.date_added}</span>
        </Tooltip>
      </div>
      <div>
        <Icon icon="ri:calendar-check-line" inline />
        <Tooltip text={`${n_days_ago(model.date_published)} days ago`}>
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
          pip install {model.pypi.split('/').pop()}
          <!-- TODO add custom CopyButton labels to remove text -->
          <CopyButton />
        </code>
      {/if}
    </section>

    <section class="links">
      <a href={model.repo} target="_blank" rel="noopener noreferrer">
        <Icon icon="octicon:mark-github" inline /> Repository
      </a>
      <a href={model.paper} target="_blank" rel="noopener noreferrer">
        <Icon icon="ion:ios-paper" inline /> Paper
      </a>
      {#if model.url}
        <a href={model.url} target="_blank" rel="noopener noreferrer">
          <Icon icon="ion:ios-globe" inline /> Documentation
        </a>
      {/if}
      <a href={model.doi} target="_blank" rel="noopener noreferrer">
        <Icon icon="academicons:doi" inline /> DOI
      </a>
      <a
        href={`${pkg.repository}/blob/-/models/${model.dirname.split('/').pop()}`}
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
        <li><strong>Type:</strong> {model.model_type}</li>
        <li><strong>Targets:</strong> {model.targets}</li>
        <li><strong>Openness:</strong> {model.openness}</li>
        <li><strong>Train Task:</strong> {model.train_task}</li>
        <li><strong>Test Task:</strong> {model.test_task}</li>
        <li>
          <strong>Trained for Benchmark:</strong>
          {model.trained_for_benchmark ? 'Yes' : 'No'}
        </li>
      </ul>
    </section>

    {#if model.training_set}
      {@const train_set_info =
        typeof model.training_set == 'string'
          ? TRAINING_SETS[model.training_set]
          : model.training_set}
      {@const { n_structures, url, title, n_materials } = train_set_info}
      {@const pretty_n_mat =
        typeof n_materials == `number` ? pretty_num(n_materials) : n_materials}
      {@const n_mat_str = n_materials ? ` from ${pretty_n_mat} materials` : ``}
      <section class="training-set">
        <h2>Training Set</h2>
        <ul>
          <li>
            <strong>Title:</strong>
            <a href={url} target="_blank" rel="noopener noreferrer">{title}</a>
          </li>
          <li>
            <strong>Number of Structures:</strong>
            <Tooltip text={n_structures.toLocaleString()}>
              <span>{pretty_num(n_structures)}</span>
            </Tooltip>
          </li>
          {#if n_materials}<li>
              <strong>Number of Materials:</strong>
              <Tooltip text={n_materials.toLocaleString()}>
                <span>{pretty_n_mat}</span>
              </Tooltip>
            </li>{/if}
        </ul>
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
      <section class="requirements">
        <h2>Package Dependencies</h2>
        <ul>
          {#each Object.entries(model.requirements) as [pkg, version]}
            {@const href = version.startsWith('http')
              ? version
              : `https://pypi.org/project/${pkg}/${version}`}
            <li>
              <strong>{pkg}</strong>:
              <a {href} target="_blank" rel="noopener noreferrer">{version}</a>
            </li>
          {/each}
        </ul>
      </section>
    {/if}

    {#if model.notes}
      {@const { description, training, missing_preds, ...otherNotes } = model.notes}
      <section class="notes">
        <h2>Notes</h2>
        {#if description}
          <h3>Description</h3>
          <p>{@html description}</p>
        {/if}
        {#if training}
          <h3>Training</h3>
          {#if typeof training === 'string'}
            <p>{@html training}</p>
          {:else}
            <ul>
              {#each Object.entries(training) as [key, value]}
                <li><strong>{key}:</strong> {value}</li>
              {/each}
            </ul>
          {/if}
        {/if}
        {#if missing_preds}
          <h3>Missing Predictions</h3>
          <p>{missing_preds}</p>
        {/if}
        {#each Object.entries(otherNotes) as [key, value]}
          <h3>{key}</h3>
          <p>{@html value}</p>
        {/each}
      </section>
    {/if}

    {#if $page.params.slug == 'mace'}
      <Mace />
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
    margin: 0;
    padding: 0;
  }

  h3 {
    margin: 1em 0;
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
</style>
