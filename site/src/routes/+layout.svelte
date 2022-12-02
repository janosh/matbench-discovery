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
</script>

<Toc headingSelector="main > :where(h2, h3, h4)" breakpoint={1250} />

{#if $page.url.pathname !== `/`}
  <a href="." aria-label="Back to index page">&laquo; back</a>
{/if}

<GitHubCorner href={repository} />

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
