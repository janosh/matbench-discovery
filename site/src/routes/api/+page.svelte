<h1 class="toc-exclude">API</h1>

{#each Object.values(import.meta.glob(`./*.md`, { eager: true })) as file}
  <svelte:component this={file?.default} />
{/each}

<style>
  /* select all but first module h1s */
  :global(h1[id^='module-']:not(:nth-of-type(2))) {
    margin-top: 2em;
  }
  :global(h1 > kbd) {
    font-size: 25pt;
  }
  :global(hr) {
    border: 1px solid gray;
    margin-bottom: 2em;
  }
</style>
