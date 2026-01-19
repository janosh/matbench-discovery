<script lang="ts">
  import { goto } from '$app/navigation'
  import { page } from '$app/state'
  import { Footer, ThemeToggle } from '$lib'
  import pkg from '$site/package.json'
  import type { Snippet } from 'svelte'
  import { CmdPalette, CopyButton, GitHubCorner, Nav } from 'svelte-multiselect'
  import { heading_anchors } from 'svelte-multiselect/heading-anchors'
  import Toc from 'svelte-toc'
  import '../app.css'

  let { children }: { children?: Snippet } = $props()

  const routes = Object.keys(import.meta.glob(`./*/+page.{svelte,md}`)).map(
    (filename) => `/` + filename.split(`/`)[1],
  )

  let url = $derived(page.url.pathname)
  let headingSelector = $derived(
    `main :is(${{ '/api': `h1, ` }[url] ?? ``}h2, h3, h4):not(.toc-exclude)`,
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
  }
  let description = $derived(
    descriptions[url ?? ``] ?? base_description,
  )
  let title = $derived(url == `/` ? `` : `${url} • `)

  const actions = Object.keys(import.meta.glob(`./**/+page.{svelte,md}`)).map(
    (filename) => {
      const parts = filename.split(`/`).filter((part) => !part.startsWith(`(`)) // remove hidden route segments
      const route = `/${parts.slice(1, -1).join(`/`)}`

      return { label: route, action: () => goto(route) }
    },
  )
</script>

<CmdPalette {actions} placeholder="Go to..." />
<CopyButton global />

<svelte:head>
  <title>{title}Matbench Discovery</title>
  <meta name="description" content={description} />
</svelte:head>

{#if ![`/`, `/models`, `/tasks/geo-opt`].includes(url)}
  <Toc
    {headingSelector}
    breakpoint={1600}
    minItems={3}
    aside_style="max-width: 22em; z-index: 1"
    nav_style="font-size: 7pt"
    --toc-mobile-width="22em"
    --toc-active-bg="transparent"
    --toc-active-color="var(--link-color)"
  />
{/if}

<GitHubCorner href={pkg.repository} />

<Nav
  routes={[`/`, ...routes.filter((route) => route != `/changelog`), [pkg.paper, `Paper`]]}
  style="left: initial; margin-block: 1em 0"
  menu_props={{ style: `gap: 1.5em; place-items: center` }}
  labels={{ '/': `Home`, '/api': `API` }}
  link_props={{
    style: `background-color: var(--nav-bg); display: inline-flex; place-items: center`,
  }}
>
  <ThemeToggle />
</Nav>

<main {@attach heading_anchors()}>
  {@render children?.()}
</main>

<Footer />

<style>
  :global(aside.toc.desktop) {
    position: fixed;
    left: calc(50vw + var(--main-max-width) / 2 + 1em);
  }
  :global(aside.toc.mobile > nav) {
    box-shadow: 0 0 20px var(--shadow);
    border: 1px solid var(--border);
    background-color: var(--page-bg);
    padding: 1ex;
    right: 1em;
    position: fixed;
    bottom: 1em;
  }
</style>
