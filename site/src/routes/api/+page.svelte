{#each Object.values(import.meta.glob(`./*.md`, { eager: true })) as file}
  <svelte:component this={file?.default} />
{/each}

<style>
  :global(hr) {
    border: none;
    margin: 3em;
  }
  :global(code) {
    line-height: 1em;
    border-radius: 4pt;
  }
</style>
