<script lang="ts">
  import { beforeUpdate } from 'svelte'
  import type { Reference } from '.'

  export let references: Reference[]
  export let ref_selector: string = `a[href^='#']`
  export let found_on_page: Reference[] = references

  function filter_refs() {
    const ref_links = document.querySelectorAll<HTMLAnchorElement>(ref_selector)
    const hrefs = Array.from(ref_links).map((ref) => ref.hash)
    found_on_page = references.filter((ref) => hrefs.includes(`#${ref.id}`))
  }
  beforeUpdate(filter_refs)
</script>

<ol>
  {#each found_on_page as { title, id, author, DOI, URL, issued }}
    <li>
      <strong {id}>{title}</strong>
      <span>
        {@html author.map((a) => `${a.given} ${a.family}`).join(`, &thinsp; `)}
      </span>
      &mdash;
      <small>
        {#if DOI}
          DOI: <a href="https://doi.org/{DOI}">{DOI}</a>
        {:else if URL}
          preprint: <a href={URL}>{URL}</a>
        {/if}
        {#if issued}
          &mdash; {issued[0].year}
        {/if}
      </small>
    </li>
  {/each}
</ol>

<style>
  ol > li {
    margin: 1ex 0;
  }
  ol > li > strong {
    display: block;
  }
</style>
