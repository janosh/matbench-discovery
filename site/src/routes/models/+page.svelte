<script lang="ts">
  import { ModelCard } from '$lib'
  import type { PageData } from './$types'

  export let data: PageData
</script>

<h1>Models</h1>

<ol>
  {#each data.models as [key, metadata], idx}
    <li>
      <ModelCard {key} data={metadata} />
    </li>
  {/each}
</ol>

<style>
  ol {
    display: grid;
    gap: 2em;
    list-style: none;
  }
  ol > li {
    background-color: rgba(255, 255, 255, 0.05);
    padding: 3pt 9pt 5pt;
  }
</style>
