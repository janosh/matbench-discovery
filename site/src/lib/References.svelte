<script lang="ts">
  import type { Reference } from '$lib'
  import type { HTMLAttributes } from 'svelte/elements'

  let {
    references,
    ref_selector = `a.ref[href^='#']`,
    found_on_page = $bindable(references),
    n_authors = 1,
    first_name_mode = `none`,
    ...rest
  }: HTMLAttributes<HTMLOListElement> & {
    references: Reference[]
    ref_selector?: string
    found_on_page?: Reference[]
    n_authors?: number
    first_name_mode?: `initial` | `full` | `none`
  } = $props()

  function filter_refs() {
    const ref_links = document.querySelectorAll<HTMLAnchorElement>(ref_selector)
    const hashes = Array.from(ref_links).map((ref) => ref.hash)
    found_on_page = references.filter((ref) => hashes.includes(`#${ref.id}`))
  }
  $effect.pre(filter_refs)
</script>

{#key found_on_page}
  <ol {...rest}>
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
