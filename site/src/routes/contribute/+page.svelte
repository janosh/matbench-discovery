<script lang="ts">
  import { data_files } from '$lib'
  import Contributing from '$root/contributing.md'
</script>

<Contributing>
  <ol slot="data-files">
    {#each Object.entries(data_files).filter(([key]) => !key.startsWith(`_`)) as [key, { url, path, html }]}
      <li style="margin-top: 1ex;">
        <strong><code>{key}</code></strong>
        (<a href={url}>{path}</a>)<br />
        {@html html}
      </li>
    {/each}
  </ol>
</Contributing>
