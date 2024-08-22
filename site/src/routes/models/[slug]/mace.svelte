<script lang="ts">
  import { browser } from '$app/environment'
  import EformParity from '$figs/mace-wbm-IS2RE-e-form-parity.svelte'
  import RawEnergyParity from '$figs/mace-wbm-IS2RE-raw-energy-parity.svelte'
</script>

<h2>MACE Analysis for WBM IS2RE</h2>

<h3>Formation Energy Parity</h3>

{#if browser}
  <EformParity style="margin: 2em 0;" />
{/if}

<h3>Raw Energy Parity</h3>

{#if browser}
  <RawEnergyParity style="margin: 2em 0;" />
{/if}
