<script lang="ts">
  import { browser } from '$app/environment'
  import RawEnergyParity from '$figs/mace-wbm-IS2RE-raw-energy-parity.svelte'
</script>

### Raw DFT Energy Parity

Total energy without [MP 2020 energy corrections](https://github.com/materialsproject/pymatgen/blob/02a4ca8aa0277b5f6db11f4de4fdbba129de70a5/pymatgen/entries/compatibility.py#L823).

{#if browser}
<RawEnergyParity style="margin: 2em 0;" />
{/if}

<style>
  p {
    text-align: center;
    margin: 0.5em auto -2em;
  }
</style>
