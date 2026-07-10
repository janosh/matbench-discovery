<script lang="ts">
  import type { PageData } from './$types'

  let { data }: { data: PageData } = $props()
</script>

{@html data.changelog?.code}
