<script lang="ts">
  import { beforeUpdate } from 'svelte'
  import type { Reference } from '.'

  export let references: Reference[]
  export let ref_selector: string = `a.ref[href^='#']`
  export let found_on_page: Reference[] = references

  function filter_refs() {
    const ref_links = document.querySelectorAll<HTMLAnchorElement>(ref_selector)
    const hashes = Array.from(ref_links).map((ref) => ref.hash)
    found_on_page = references.filter((ref) => hashes.includes(`#${ref.id}`))
  }
  beforeUpdate(filter_refs)
</script>

{#key found_on_page}
  <ol>
    {#each found_on_page as { title, id, author, DOI, URL: href, issued } (id)}
      <li>
        <strong {id}>{title}</strong>
        <span>
          {@html author.map((a) => `${a.given} ${a.family}`).join(`, &thinsp; `)}
        </span>
        &mdash;
        <small>
          {#if DOI}
            <a href="https://doi.org/{DOI}">{DOI}</a>
          {:else if href}
            <a {href}>{href}</a>
          {/if}
          {#if issued}
            &mdash; {issued[0].year}
          {/if}
        </small>
      </li>
    {/each}
  </ol>
{/key}

<style>
  ol {
    padding: 0 0 0 1em;
  }
  ol > li {
    margin: 1ex 0;
  }
  ol > li > strong {
    display: block;
  }
  ol > li > :is(small, span) {
    font-weight: lighter;
  }
</style>
