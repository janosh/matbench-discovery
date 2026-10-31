<script lang="ts">
  import type { Component } from 'svelte'
  import Intro from './intro.md'

  const modules = Object.entries(
    import.meta.glob<Component>(`./modules/*.md`, {
      eager: true,
      import: 'default',
    }),
  ).toSorted(([path_a], [path_b]) => path_a.localeCompare(path_b))
</script>

<div class="api">
  <h1 id="api" class="toc-exclude">API</h1>

  <Intro />

  {#each modules as [path, MdFile] (path)}
    <MdFile />
  {/each}
</div>

<style>
  /* select all but first module h1s */
  .api :global(h1[id^='module-']:not(:nth-of-type(2))) {
    margin-top: 2em;
  }
  .api :global(h1 > kbd) {
    font-size: 25pt;
  }
  .api :global(hr) {
    border: 1px solid var(--border);
    margin-bottom: 2em;
  }
</style>
