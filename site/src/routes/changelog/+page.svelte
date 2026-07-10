<script lang="ts">
  let { data }: { data: { changelog?: { code?: string } } } = $props()
</script>

{@html data.changelog?.code}
