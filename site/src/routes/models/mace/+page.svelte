<script lang="ts">
  import { browser } from '$app/environment'
  import EformParity from '$figs/mace-wbm-IS2RE-e-form-parity.svelte'
  import RawEnergyParity from '$figs/mace-wbm-IS2RE-raw-energy-parity.svelte'
</script>

<h1>MACE Additional Analysis for WBM IS2RE</h1>

<h2>Formation Energy Parity</h2>

{#if browser}
  <EformParity style="margin: 2em 0;" />
{/if}

<h2>Raw Energy Parity</h2>

{#if browser}
  <RawEnergyParity style="margin: 2em 0;" />
{/if}
