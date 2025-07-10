<script lang="ts">
  import { goto } from '$app/navigation'
  import { page } from '$app/state'
  import { Footer, Nav } from '$lib'
  import pkg from '$site/package.json'
  import { type Snippet } from 'svelte'
  import { CmdPalette, CopyButton, GitHubCorner } from 'svelte-multiselect'
  import Toc from 'svelte-toc'
  import '../app.css'

  interface Props {
    children?: Snippet
  }
  let { children }: Props = $props()

  const routes = Object.keys(import.meta.glob(`./*/+page.{svelte,md}`)).map(
    (filename) => `/` + filename.split(`/`)[1],
  )

  let url = $derived(page.url.pathname)
  let headingSelector = $derived(
    `main :is(${{ '/api': `h1, ` }[url] ?? ``}h2, h3, h4):not(.toc-exclude)`,
  )

  let description = $derived(
    {
      '/': `Benchmarking machine learning energy models for materials discovery.`,
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
    }[url ?? ``],
  )
  $effect(() => {
    if (url && !description) console.warn(`No description for url=${url}`)
  })
  let title = $derived(url == `/` ? `` : `${url} • `)

  const actions = Object.keys(import.meta.glob(`./**/+page.{svelte,md}`)).map(
    (filename) => {
      const parts = filename.split(`/`).filter((part) => !part.startsWith(`(`)) // remove hidden route segments
      const route = `/${parts.slice(1, -1).join(`/`)}`

      return { label: route, action: () => goto(route) }
    },
  )
  $effect.pre(() => {
    if (page.url.pathname == `/models`) {
      document.documentElement.style.setProperty(`--main-max-width`, `90em`)
    } else {
      document.documentElement.style.setProperty(`--main-max-width`, `50em`)
    }
  })
</script>

<CmdPalette {actions} placeholder="Go to..." />
<CopyButton global />

<svelte:head>
  <title>{title}Matbench Discovery</title>
  <meta name="description" content={description} />
</svelte:head>

{#if ![`/`, `/models`].includes(url)}
  <Toc
    {headingSelector}
    breakpoint={1600}
    minItems={3}
    aside_style="min-width: 25em; left: calc(50vw + var(--main-max-width) / 2 + 160px); line-height: 1.8;"
    open_button_style="right: 2em;"
    --toc-mobile-bg="rgba(0, 0, 0, 0.3)"
    --toc-mobile-width="22em"
    --toc-active-bg="none"
  />
{/if}

<GitHubCorner href={pkg.repository} />

<Nav
  routes={[
    [`/home`, `/`],
    ...routes.filter((route) => route != `/changelog`),
    [`/paper`, pkg.paper],
  ]}
  style="padding: 0 var(--main-padding)"
/>

<main>
  {@render children?.()}
</main>

<Footer />
