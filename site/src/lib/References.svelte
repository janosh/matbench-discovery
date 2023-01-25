<script lang="ts">
  import type { Reference } from '.'

  export let references: Reference[]
</script>

<ol>
  {#each references as { title, id, author, DOI, URL, issued }}
    <li>
      <strong {id}>{title}</strong>
      <p>
        {@html author.map((a) => `${a.given} ${a.family}`).join(`, &thinsp; `)}
      </p>
      <p>
        {#if DOI}
          DOI: <a href="https://doi.org/{DOI}">{DOI}</a>
        {:else if URL}
          preprint: <a href={URL}>{URL}</a>
        {/if}
        {#if issued}
          - {issued[0].year}
        {/if}
      </p>
    </li>
  {/each}
</ol>

<style>
  ol > li {
    margin: 1ex 0;
  }
  ol > li > p {
    margin: 0;
  }
</style>
