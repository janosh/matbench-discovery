<h1 class="toc-exclude">API</h1>

{#each Object.entries(import.meta.glob(`./*.md`, { eager: true })) as
  [path, { default: MdFile }]
  (path)
}
  <MdFile />
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
