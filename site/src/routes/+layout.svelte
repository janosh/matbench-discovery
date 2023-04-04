<script lang="ts">
  import { goto } from '$app/navigation'
  import { page } from '$app/stores'
  import { Footer, Nav } from '$lib'
  import { repository } from '$site/package.json'
  import { CmdPalette } from 'svelte-multiselect'
  import Toc from 'svelte-toc'
  import { GitHubCorner, PrevNext } from 'svelte-zoo'
  import '../app.css'

  const routes = Object.keys(import.meta.glob(`./*/+page.{svx,svelte,md}`)).map(
    (filename) => `/` + filename.split(`/`)[1]
  )

  $: url = $page.url.pathname
  $: headingSelector = `main :is(${
    url === `/api` ? `h1, ` : ``
  }h2, h3, h4):not(.toc-exclude)`

  $: description = {
    '/': `Benchmarking machine learning energy models for materials discovery.`,
    '/about-the-data': `Details about provenance, chemistry and energies in the benchmark's train and test set.`,
    '/about-the-data/tmi': `Too much information on the benchmark's data.`,
    '/api': `API docs for the Matbench Discovery PyPI package.`,
    '/contribute': `Steps for contributing a new model to the benchmark.`,
    '/models': `Details on each model sortable by metrics.`,
    '/paper': `The paper released with the Matbench Discovery benchmark.`,
    '/paper/iclr-ml4mat': `Extended abstract submitted to the ICLR ML4Materials workshop.`,
    '/si': `Supplementary information including interesting but non-essential plots.`,
  }[url ?? ``]
  if (url && !description) console.warn(`No description for url=${url}`)
  $: title = url == `/` ? `` : `${url} • `

  const actions = Object.keys(import.meta.glob(`./**/+page.{svx,svelte,md}`)).map(
    (filename) => {
      const parts = filename.split(`/`).filter((part) => !part.startsWith(`(`)) // remove hidden route segments
      const route = `/${parts.slice(1, -1).join(`/`)}`

      return { label: route, action: () => goto(route) }
    }
  )
</script>

<CmdPalette {actions} placeholder="Go to..." />

<svelte:head>
  <title>{title}Matbench Discovery</title>
  <meta name="description" content={description} />
</svelte:head>

<Toc {headingSelector} breakpoint={1250} minItems={3} />

{#if url !== `/`}
  <a href="/" aria-label="Back to index page">&laquo; home</a>
{/if}

<GitHubCorner href={repository} />

<main>
  <Nav {routes} />

  <slot />

  <PrevNext
    items={routes}
    current="/{url?.split(`/`)[1]}"
    let:item={href}
    style="margin-top: 4em;"
  >
    <a {href} class="link" slot="next">{href} &raquo;</a>
    <a {href} class="link" slot="prev">&laquo; {href}</a>
  </PrevNext>
</main>

<Footer />

<style>
  main {
    margin: auto;
    margin-bottom: 3em;
    width: 100%;
    max-width: 50em;
  }
  a[href='/'] {
    font-size: 14pt;
    background-color: rgba(255, 255, 255, 0.1);
    padding: 1pt 5pt;
    border-radius: 3pt;
    transition: 0.2s;
    width: max-content;
    box-sizing: border-box;
    margin: 1em 0 0 1em;
  }
  @media (min-width: 900px) {
    a[href='/'] {
      position: absolute;
      top: 1em;
      left: 1em;
    }
  }
  a[href='/']:hover {
    background-color: rgba(255, 255, 255, 0.2);
  }
</style>
