<script>
  import data_files from '$pkg/data-files.yml'
</script>

## Downloading data

Public benchmark inputs and reference files are available from [Figshare](https://figshare.com/articles/dataset/22715158). The `DataFiles` registry downloads and caches individual files locally:

```py
from matbench_discovery.enums import DataFiles

wbm_structures = DataFiles.wbm_initial_atoms.path
phonon_structures = DataFiles.phonondb_pbe_103_structures.path
md_trajectories = DataFiles.dynamat_v1_0_md_trajectories.path
diatomic_curves = DataFiles.diatomics_dft_reference.path
```

Each task’s test-set section links its public inputs and explains reference availability.

### All public data files

<ol class="data-files-list">
{#each Object.entries(data_files).filter(([key]) => !key.startsWith(`_`)) as [key, { url, path, html }]}
    <li style="margin-top: 1ex;">
    <strong><code>{key}</code></strong>
    {#if url}
      (<a href={url}>{path}</a>)
    {:else}
      (<code>{path}</code>, not yet published)
    {/if}<br />
    {@html html}
    </li>
{/each}
</ol>
