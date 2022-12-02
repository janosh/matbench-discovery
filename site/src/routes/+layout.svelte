<script lang="ts">
  import { page } from '$app/stores'
  import { repository } from '$site/package.json'
  import GitHubCorner from 'svelte-github-corner'
  import Toc from 'svelte-toc'
  import '../app.css'

  $: is_current = (path: string) => {
    if (path === $page.url.pathname) return `page`
    if (path !== `/` && $page.url.pathname.includes(path)) return `page`
    return undefined
  }

  const routes = Object.keys(import.meta.glob(`./*/+page.{svx,svelte,md}`)).map(
    (filename) => filename.split(`/`)[1]
  )
</script>

<Toc headingSelector="main > :is(h2, h3, h4):not(.toc-exclude)" breakpoint={1250} />

{#if $page.url.pathname !== `/`}
  <a href="." aria-label="Back to index page">&laquo; back</a>
{/if}

<GitHubCorner href={repository} />

<nav>
  {#each routes as route, idx}
    {#if idx > 0}<strong>&bull;</strong>{/if}
    <a href={route} aria-current={is_current(route)}>{route}</a>
  {/each}
</nav>

<main>
  <slot />
</main>

<style>
  :global(aside.toc.desktop) {
    position: fixed;
    top: 3em;
    left: calc(50vw + 50em / 2);
    max-width: 16em;
  }
  nav {
    display: flex;
    gap: 1em 1ex;
    place-content: center;
    margin: 1em auto 3em;
    max-width: 45em;
    flex-wrap: wrap;
    font-size: 1.1em;
  }
  nav > a {
    padding: 0 4pt;
    background-color: rgba(255, 255, 255, 0.1);
    border-radius: 3pt;
    transition: 0.2s;
  }
  nav > a[aria-current='page'] {
    color: mediumseagreen;
  }
  a[href='.'] {
    font-size: 15pt;
    position: absolute;
    top: 2em;
    left: 2em;
    background-color: rgba(255, 255, 255, 0.1);
    padding: 1pt 5pt;
    border-radius: 3pt;
    transition: 0.2s;
  }
  a[href='.']:hover {
    background-color: rgba(255, 255, 255, 0.2);
  }
</style>
