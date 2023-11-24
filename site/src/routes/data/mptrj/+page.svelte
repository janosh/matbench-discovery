<script lang="ts">
  import figshare_urls from '$root/data/figshare/1.0.0.json'

  const plots = import.meta.glob(`$figs/mp-trj-*.svelte`, {
    eager: true,
    import: `default`,
  })

  const title_map: Record<string, string> = {
    'e-form': `Formation Energy`,
    magmoms: `Magnetic Moments`,
  }
</script>

<h1>MPtrj Target Distributions</h1>

The MPtrj universal potential training set
<a href={figshare_urls.mptrj.article}>available on figshare</a>
contains 1,580,395 structures, 1,580,395 energies, 7,944,833 magnetic moments, 49,295,660 forces,
and 14,223,555 stresses that were used to train the open-source interatomic potentials tested
on Matbench Discovery.

<ul>
  {#each Object.entries(plots) as [name, figure]}
    {@const title = name.split(`mp-trj-`)[1].split(`-hist.svelte`)[0]}
    <h2>{title_map[title] ?? title}</h2>
    <svelte:component this={figure} style="width: 100%;" {title} />
  {/each}
</ul>

<style>
  ul {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(400px, 1fr));
    gap: 1em;
  }
  ul h2 {
    text-transform: capitalize;
    margin: 1em auto 0;
  }
</style>
