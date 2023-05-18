<script lang="ts">
  import Select from 'svelte-multiselect'

  export let disabled: boolean = false

  const figs = import.meta.glob(
    `$figs/scatter-largest-errors-models-mean-vs-true-hull-dist-*.svelte`,
    { eager: true }
  )

  let selected: string[] = [Object.keys(figs)[0]]
</script>

<Select options={Object.keys(figs)} bind:selected minSelect={1} maxSelect={1} {disabled}>
  <span let:option slot="selected">
    {option.split(`-`).slice(-1)[0].split(`.`)[0]}
  </span>
  <span let:option slot="option">
    {option.split(`-`).slice(-1)[0].split(`.`)[0]}
  </span>
</Select>

<svelte:component this={figs[selected[0]]?.default} />

<style>
  span {
    text-transform: capitalize;
  }
</style>
