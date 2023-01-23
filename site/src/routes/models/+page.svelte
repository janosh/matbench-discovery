<script lang="ts">
  import { ModelCard } from '$lib'
  import type { PageData } from './$types'

  export let data: PageData
</script>

<h1 class="pull-left">Models</h1>

<ol class="pull-left">
  {#each data.models as [key, metadata]}
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
    grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
  }
  ol > li {
    background-color: rgba(255, 255, 255, 0.05);
    padding: 3pt 10pt 7pt;
    border-radius: 3pt;
  }
</style>
