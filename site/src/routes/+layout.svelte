<script lang="ts">
  import { goto } from '$app/navigation'
  import { page } from '$app/state'
  import { Footer } from '$lib'
  import { MODELS } from '$lib/models.svelte'
  import pkg from '$site/package.json'
  import type { Snippet } from 'svelte'
  import { CmdPalette, CopyButton, GitHubCorner, Nav, ThemeToggle } from 'svelte-multiselect'
  import { heading_anchors } from 'svelte-multiselect/heading-anchors'
  import Toc from 'svelte-toc'
  // oxlint-disable-next-line no-unassigned-import
  import '../app.css'

  let { children }: { children?: Snippet } = $props()
  let toc_desktop = $state(true)

  const routes = Object.keys(import.meta.glob(`./*/+page.{svelte,md}`)).map(
    (filename) => `/${filename.split(`/`)[1]}`,
  )

  let url = $derived(page.url.pathname)
  let headingSelector = $derived(
    `main :is(${url === `/api` ? `h1, ` : ``}h2, h3, h4)`,
  )

  const base_description =
    `Matbench Discovery - Benchmarking machine learning energy models for materials discovery.`
  const descriptions: Record<string, string> = {
    '/': base_description,
    '/data':
      `Details about provenance, chemistry and energies in the benchmark's train and test set.`,
    '/data/tmi': `Too much information on the benchmark's data.`,
    '/api': `API docs for the Matbench Discovery PyPI package.`,
    '/contribute': `Steps for contributing a new model to the benchmark.`,
    '/models': `Details on each model sortable by metrics.`,
    '/tasks/diatomics': `Metrics and analysis of predicting diatomic energies.`,
    '/tasks/phonons':
      `Metrics and analysis of predicting phonon modes and frequencies.`,
    '/tasks/geo-opt': `Metrics and analysis of predicting ground state geometries.`,
    '/tasks/md': `Metrics and analysis of molecular dynamics observables vs ab-initio reference trajectories.`,
  }
  let description = $derived(
    descriptions[url ?? ``] ?? base_description,
  )
  let title = $derived(url === `/` ? `` : `${url} • `)

  const actions = Object.keys(import.meta.glob(`./**/+page.{svelte,md}`))
    .filter((filename) => !filename.includes(`[`))
    .map((filename) => {
      const parts = filename.split(`/`).filter((part) => !part.startsWith(`(`)) // Remove hidden route segments
      const route = `/${parts.slice(1, -1).join(`/`)}`

      return { label: route, action: () => goto(route) }
    })
    .concat(
      MODELS.map((model) => ({
        label: `/models/${model.model_key}`,
        action: () => goto(`/models/${model.model_key}`),
      }))
    )
</script>

<CmdPalette
  {actions}
  placeholder="Go to..."
  dialog_style="top: 15vh; bottom: auto; overflow: visible"
/>
<CopyButton global />

<svelte:head>
  <title>{title}Matbench Discovery</title>
  <meta name="description" content={description} />
</svelte:head>

{#if ![`/`, `/models`, `/tasks/geo-opt`].includes(url)}
  <Toc
    {headingSelector}
    breakpoint={1350}
    minItems={3}
    hideOnIntersect="section.full-bleed .table-container, .bleed-1400"
    bind:desktop={toc_desktop}
    asideProps={{
      style: toc_desktop ? `max-width: 22em; position: fixed; left: calc(50vw + var(--main-max-width) / 2); top: 8em;` : `z-index: 1;`,
    }}
    navProps={{
      style: toc_desktop ? `font-size: 7pt;` : `font-size: 7pt; z-index: 10; padding: 1em;`,
    }}
    titleProps={{ style: `margin: 3pt` }}
    --toc-active-color="var(--link-color)"
    --toc-padding="1em 1em 0 1.5em"
    --toc-mobile-width="min(80vw, 30em)"
    --toc-mobile-border="1px solid var(--border)"
    --toc-mobile-shadow="0 0 20px var(--shadow)"
  />
{/if}

<GitHubCorner href={pkg.repository} />

<Nav
  {page}
  routes={[`/`, ...routes.filter((route) => route != `/changelog`), [pkg.paper, `Paper`]]}
  style="margin-block: 1em 0"
  menu_props={{ style: `gap: 1.5em` }}
  labels={{ '/': `Home`, '/api': `API` }}
  link_props={{
    style: `padding: 0 3pt`,
  }}
  --nav-link-active-color="var(--link-color)"
>
  <ThemeToggle />
</Nav>

<main {@attach heading_anchors()}>
  {@render children?.()}
</main>

<Footer />
