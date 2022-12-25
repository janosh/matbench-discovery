<script lang="ts">
  import Nav from '$lib/Nav.svelte'

  const routes = Object.keys(import.meta.glob(`./*.md`)).map(
    (filename) => filename.split(`/`)[1].split(`.`)[0]
  )
</script>

<Nav {routes} />

<slot />
