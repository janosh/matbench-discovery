<script lang="ts">
  import { References } from '$lib'
  import cite from '$root/citation.cff'
  import { pretty_num } from 'elementari'
  import { references } from './references.yaml'

  export let data
</script>

<h1>{cite.title}<br /><small>{cite.subtitle}</small></h1>

<address>
  <span>
    {#each cite.authors as auth, idx}
      {#if idx > 0},
      {/if}
      <a href={auth.github ?? `https://orcid.org/${auth.orcid}`}>
        {auth[`given-names`]} {auth[`family-names`]}</a
      ><sup>{`${auth.affil_key}`.replaceAll(` `, ``)}</sup>
    {/each}
  </span>
  <span style="font-weight: lighter; font-size: 0.96em;">
    {@html cite.affiliations.map((affil, idx) => `${idx + 1}. ${affil}`).join(`<br/>`)}
  </span>
  <span style="font-weight: lighter;">{cite[`date-released`]}</span>
</address>

<div>
  <slot />
</div>

<h2>References</h2>

<References {references} />

<small style="float: right;">
  <code>{pretty_num(data.word_count)}</code> words (<code>
    ~{Math.floor(data.word_count / 150)}</code
  > min)
</small>

<style>
  address,
  h1 > small {
    text-align: center;
    text-wrap: balance;
    font-style: normal;
  }
  address sup {
    white-space: nowrap;
    padding: 0 0 0 1px;
  }
  h1 > small {
    display: block;
    font-size: 0.7em;
  }
  address span {
    margin: 1em;
    display: block;
  }
  div :global(summary.abstract) {
    font-weight: 300;
    font-size: 0.95em;
  }
  div :global(summary.abstract)::before {
    content: 'Abstract';
    font-weight: bold;
    font-size: larger;
  }
  /* references */
  div :global(.ref) {
    text-transform: capitalize;
  }
  /* auto-number HTML headings */
  div :global(h2) {
    counter-increment: h2;
    counter-reset: h3 h4;
  }
  div :global(h2::before) {
    content: counter(h2);
    margin-right: 10pt;
  }
  div :global(h3) {
    counter-increment: h3;
    counter-reset: h4;
  }
  div :global(h3::before) {
    content: counter(h2) '.' counter(h3);
    margin-right: 10pt;
  }
</style>
