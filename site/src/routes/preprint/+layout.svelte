<script lang="ts">
  import { References } from '$lib'
  import cite from '$root/citation.cff'
  import { pretty_num } from 'elementari'
  import { references } from './references.yaml'

  export let data

  const authors = cite.authors.map(
    (auth) => `${auth[`given-names`]} ${auth[`family-names`]}<sup>${auth.affil_key}</sup>`
  )
</script>

<h1>{cite.title}<br /><small>{cite.subtitle}</small></h1>

<address>
  <span>
    {@html authors.join(`, `)}
  </span>
  <span>
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
    ~{Math.floor(data.word_count / 200)}</code
  > min)
</small>

<style>
  address {
    text-align: center;
  }
  address span {
    margin: 1em;
    display: block;
  }
  div :global(summary.abstract) {
    font-weight: 300;
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
