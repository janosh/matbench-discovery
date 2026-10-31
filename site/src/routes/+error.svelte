<script lang="ts">
  import { page } from '$app/state'
  import pkg from '#site/package.json'
  import { ErrorPage } from 'svelte-widgets'
</script>

<ErrorPage
  status={page.status}
  message={page.error?.message}
  name={pkg.name}
  bugs_url="{pkg.repository}/issues"
  status_page_url="https://githubstatus.com"
  logo="/favicon.svg"
/>
