<script lang="ts">
  import Select from 'svelte-multiselect'

  export let disabled: boolean = false

  const figs = import.meta.glob(
    `$figs/scatter-largest-errors-models-mean-vs-true-hull-dist-*.svelte`,
    { eager: true, import: 'default' }
  )

  let selected: string[] = [Object.keys(figs)[0]]

  const get_chemistry = (str: string) => str.split(`-`).slice(-1)[0].split(`.`)[0]
</script>

<Select
  options={[...Object.keys(figs), `Grid`]}
  bind:selected
  minSelect={1}
  maxSelect={1}
  {disabled}
>
  <span let:option slot="selected">
    {get_chemistry(option)}
  </span>
  <span let:option slot="option">
    {get_chemistry(option)}
  </span>
</Select>

{#if selected[0] == `Grid`}
  <ul>
    {#each Object.keys(figs) as key}
      {@const chem = get_chemistry(key)}
      <li>
        <h3>{chem}</h3>
        <svelte:component this={figs[key]} />
      </li>
    {/each}
  </ul>
{:else}
  <svelte:component this={figs[selected[0]]} />
{/if}

<style>
  span {
    text-transform: capitalize;
  }
  ul {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(400px, 1fr));
    padding: 0;
    list-style: none;
    margin-left: calc(0.6 * (-50vw + 50cqw));
  }
  /* hide plotly titles */
  ul :global(g.g-gtitle) {
    display: none;
  }
  h3 {
    margin: 1em 0 -3em;
    text-align: center;
    text-transform: capitalize;
  }
</style>
