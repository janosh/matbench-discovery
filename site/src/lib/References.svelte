<script lang="ts">
  import type { Reference } from '$lib'
  import { beforeUpdate } from 'svelte'

  export let references: Reference[]
  export let ref_selector: string = `a.ref[href^='#']`
  export let found_on_page: Reference[] = references
  export let n_authors: number = 1
  export let first_name_mode: `initial` | `full` | `none` = `none`

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
        <p {id}>{title}</p>
        <span>
          {@html author
            .slice(0, n_authors)
            .map(({ given, family }) => {
              const first_name = {
                initial: `${given[0]}. `,
                full: `${given} `,
                none: ``,
              }[first_name_mode]
              return `${first_name ?? ``}${family}`
            })
            .join(`,&thinsp; `)}
          {#if author.length > n_authors}
            <em>et al.</em>
          {/if}
        </span>
        <small>
          &mdash;
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
  ol > li > p {
    margin: 0;
  }
  ol > li > :is(small, span) {
    font-weight: lighter;
  }
</style>
