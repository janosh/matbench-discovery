<script lang="ts">
  import type { PageData } from './$types'
  import Intro from './intro.md'

  let { data }: { data: PageData } = $props()
</script>

<div class="api">
  <h1 id="api" class="toc-exclude">API</h1>

  <Intro />

  {#each data.modules as html, idx (idx)}
    {@html html}
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
