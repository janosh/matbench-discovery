<script lang="ts">
  import { goto } from '$app/navigation'
  import { onMount } from 'svelte'

  // Static hosting cannot see URL fragments; resolve old dataset anchors in the browser.
  const destinations = new Map([
    [`wbm`, `/benchmarks/discovery#test-set`],
    [`phonondb`, `/benchmarks/phonons#test-set`],
    [`dynamat`, `/benchmarks/md#test-set`],
    [`diatomics`, `/benchmarks/diatomics#test-set`],
    [`training-data`, `/data/sets`],
  ])
  onMount(() => {
    const anchor = location.hash.slice(1).replace(/-title$/, ``)
    const target = new URL(
      destinations.get(anchor) ?? `/benchmarks${location.hash}`,
      location.origin,
    )
    target.search = location.search
    void goto(`${target.pathname}${target.search}${target.hash}`, { replaceState: true })
  })
</script>

<svelte:head>
  <link
    rel="canonical"
    href="https://matbench-discovery.materialsproject.org/benchmarks"
  />
</svelte:head>
<p>
  Benchmark data is now part of <a href="/benchmarks">Benchmarks</a>. Each task includes
  its test set, reference data, and downloads.
</p>
