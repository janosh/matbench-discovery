<script lang="ts">
  import { goto } from '$app/navigation'
  import { page } from '$app/state'
  import DATASETS from '$data/datasets.yml'
  import { title_case } from '$lib/labels'
  import { bind_score_weights, MODELS } from '$lib/models.svelte'
  import { bind_comparison_url, comparison } from '$lib/model-comparison.svelte'
  import { get_error_message } from '$lib/asset-loader'
  import {
    CommandMenu,
    CopyButton,
    FindBar,
    Footer,
    GitHubCorner,
    Icon,
    Nav,
    ThemeToggle,
    Toc,
  } from 'svelte-widgets'
  import type { CmdAction, FooterLink } from 'svelte-widgets'
  import { Changelog, Email, GitHub, RSS, Search, TextSearch } from 'svelte-widgets/icons'
  import MODELING_TASKS from '$pkg/modeling-tasks.yml'
  import pkg from '$site/package.json'
  import { tick, type Snippet } from 'svelte'
  import { heading_anchors } from 'svelte-widgets/heading-anchors'
  // oxlint-disable-next-line no-unassigned-import
  import '../app.css'

  let { children }: { children?: Snippet } = $props()
  let toc_desktop = $state(true)
  let find_open = $state(false)
  let search_open = $state(false)
  let find_bar = $state<ReturnType<typeof FindBar>>()
  let main_element = $state<HTMLElement>()
  bind_comparison_url()
  bind_score_weights()
  // Load on first use, then retain the dialog's axis/picker state across reopenings.
  let comparison_module =
    $state.raw<Promise<typeof import('$lib/model/ModelComparison.svelte')>>()
  $effect(() => {
    if (comparison.open) comparison_module ??= import(`$lib/model/ModelComparison.svelte`)
  })

  const footer_links: FooterLink[] = [
    { href: `${pkg.repository}/issues`, label: `Issues`, icon: GitHub },
    {
      href: `mailto:janosh.riebesell@gmail.com?subject=Matbench Discovery`,
      label: `Contact`,
      icon: Email,
    },
    { href: `/changelog`, label: `Changelog`, icon: Changelog },
    {
      href: `/rss.xml`,
      label: `RSS`,
      icon: RSS,
      title: `Be notified of new model submissions`,
    },
  ]

  // show full task titles from modeling-tasks.yml instead of capitalized URL slugs
  const task_links = Object.entries(MODELING_TASKS)
    .filter(([key]) => key !== `cps`)
    .map(([key, task]) => ({
      href: `/benchmarks/${key.replaceAll(`_`, `-`)}`,
      label: task.label,
    }))
    .toSorted((first_link, second_link) =>
      first_link.href.localeCompare(second_link.href),
    )
  const nav_links = [
    { href: `/`, label: `Home` },
    { href: `/benchmarks`, label: `Benchmarks`, children: task_links },
    { href: `/models`, label: `Models` },
    { href: `/api`, label: `API` },
    { href: `/contribute`, label: `Contribute` },
    { href: `/data/sets`, label: `Datasets` },
    { href: pkg.paper, label: `Paper` },
  ]

  let url = $derived(page.url.pathname)
  let heading_selector = $derived(`main :is(${url === `/api` ? `h1, ` : ``}h2, h3, h4)`)
  let find_enabled = $derived([`/api`, `/changelog`, `/contribute`].includes(url))

  const open_find = async (): Promise<void> => {
    find_open = true
    await tick()
    find_bar?.focus_input()
  }

  // Reset so leaving and returning to a searchable route does not reopen the bar.
  $effect(() => {
    if (!find_enabled) find_open = false
  })

  const base_description = `Matbench Discovery - ${pkg.description}`
  const descriptions: Record<string, string> = {
    '/': base_description,
    '/data': `Benchmark test sets have moved to the benchmark task pages.`,
    '/data/sets': `Explore datasets, access terms, licenses, and model training data.`,
    '/data/tmi': `Additional chemical-diversity analysis of the WBM test set.`,
    '/api': `API docs for the Matbench Discovery PyPI package.`,
    '/contribute': `Steps for contributing a new model to the benchmark.`,
    '/models': `Details on each model sortable by metrics.`,
    '/benchmarks': `Benchmark tasks, test sets, reference data, and model results for ML force fields.`,
    '/benchmarks/discovery': `Metrics and analysis of crystal stability prediction on the WBM test set.`,
    '/benchmarks/discovery/tmi': `Detailed diagnostics for the crystal discovery task.`,
    '/benchmarks/diatomics': `Metrics and analysis of predicting diatomic energies.`,
    '/benchmarks/phonons': `Metrics and analysis of predicting phonon modes and frequencies.`,
    '/benchmarks/geo-opt': `Metrics and analysis of predicting ground state geometries.`,
    '/benchmarks/md': `Metrics and analysis of molecular dynamics observables vs ab-initio reference trajectories.`,
  }
  let description = $derived(descriptions[url] ?? base_description)
  let title = $derived(url === `/` ? `` : `${url} • `)

  const page_labels = Object.fromEntries(
    [
      ...nav_links,
      ...task_links,
      { href: `/data/tmi`, label: `WBM chemical diversity` },
      { href: `/benchmarks/discovery/tmi`, label: `Discovery diagnostics` },
      { href: `/benchmarks/diatomics/tmi`, label: `Diatomic DFT reference spin states` },
    ].map(({ href, label }) => [href, label]),
  )
  const page_routes = Object.keys(import.meta.glob(`./**/+page.{svelte,md}`))
    .filter((filename) => !filename.includes(`[`))
    .map((filename) => {
      const parts = filename.split(`/`).filter((part) => !part.startsWith(`(`))
      return `/${parts.slice(1, -1).join(`/`)}`
    })
    .filter((route) => route !== `/data`)
  const link_action = (
    href: string,
    label: string,
    group: string,
    ...keywords: string[]
  ): CmdAction => ({
    id: href,
    label,
    group,
    keywords: [href, ...keywords],
    action: () => goto(href),
  })
  const actions: CmdAction[] = [
    ...page_routes.map((route) =>
      link_action(
        route,
        page_labels[route] ??
          route.split(`/`).filter(Boolean).map(title_case).join(` · `),
        `Pages`,
      ),
    ),
    ...MODELS.map(({ model_key, model_name }) =>
      link_action(`/models/${model_key}`, model_name, `Models`),
    ),
    ...Object.entries(DATASETS).map(([key, { name, slug }]) =>
      link_action(`/data/${slug}`, name, `Datasets`, key),
    ),
  ]
</script>

<CommandMenu
  {actions}
  bind:open={search_open}
  placeholder="Search models, datasets, and pages…"
  input_style="padding: 0.3em 0.6em"
  style="overflow: hidden"
  li_option_style="padding: 0.15em 0.6em"
  li_group_header_style="padding: 0.2em 0.6em"
  dialog_props={{ style: `top: 15vh` }}
/>
<CopyButton global />

<svelte:head>
  <title>{title}Matbench Discovery</title>
  <meta name="description" content={description} />
</svelte:head>

{#if ![`/`, `/models`, `/benchmarks/geo-opt`].includes(url)}
  <Toc
    {heading_selector}
    dynamic
    breakpoint={1350}
    min_items={3}
    hide_on_intersect="section.full-bleed .table-container, .bleed-1400"
    bind:desktop={toc_desktop}
    aside_props={{
      style: toc_desktop
        ? `max-width: 22em; position: fixed; left: calc(50vw + var(--main-max-width) / 2); top: 8em;`
        : `z-index: 1;`,
    }}
    nav_props={{
      style: toc_desktop
        ? `font-size: 7pt;`
        : `font-size: 7pt; z-index: 10; padding: 1em;`,
    }}
    title_props={{ style: `margin: 3pt` }}
    --toc-title-font-weight="600"
    --toc-li-color="var(--text-color)"
    --toc-active-color="var(--link-color)"
    --toc-padding="1em 1em 0 1.5em"
    --toc-mobile-width="min(80vw, 30em)"
    --toc-mobile-border="1px solid var(--border)"
    --toc-mobile-shadow="0 0 20px var(--shadow)"
  />
{/if}

<GitHubCorner href={pkg.repository} id="github-corner" />

<Nav
  pathname={url}
  routes={nav_links}
  style="margin-block: 1em 0"
  --nav-item-padding="0 3pt"
  --nav-dropdown-link-padding="2pt 4pt"
  --nav-link-active-color="var(--link-color)"
  --nav-mobile-z-index="50"
  --nav-toggle-btn-z-index="50"
>
  <button
    aria-label="Search models, datasets, and pages"
    class="find-page"
    onclick={() => (search_open = true)}
    title="Search (⌘K / Ctrl+K)"
    type="button"
  >
    <Icon icon={Search} />
  </button>
  {#if find_enabled}
    <button
      aria-label="Find in page"
      class="find-page"
      onclick={open_find}
      title="Find in page"
      type="button"
    >
      <Icon icon={TextSearch} />
    </button>
  {/if}
  <ThemeToggle />
</Nav>

<main
  bind:this={main_element}
  {@attach heading_anchors({ selector: `h1, h2, h3, h4, h5, h6` })}
>
  {#if find_open}
    <FindBar
      bind:this={find_bar}
      root={main_element}
      on_close={() => (find_open = false)}
      also_ignore="[role='search']"
      style="position: sticky; inset-inline: auto 0; margin-inline-start: auto; top: 0.5rem;"
    />
  {/if}
  {@render children?.()}
</main>

{#await comparison_module then module}
  {#if module}<module.default />{/if}
{:catch error}
  {#if comparison.open}
    <div style="padding: 1rem; text-align: center">
      <p role="alert">Could not load model comparison.</p>
      <button onclick={() => location.reload()}>Reload comparison</button>
      <button onclick={() => (comparison.open = false)}>Dismiss</button>
      <details><summary>Asset details</summary>{get_error_message(error)}</details>
    </div>
  {/if}
{/await}

<Footer links={footer_links} style="--footer-bg: var(--nav-bg)">
  <img src="/favicon.svg" alt="Logo" width="30px" style="vertical-align: middle" />
  &ensp;{pkg.title} &ensp; | &ensp; ©
  <a href={pkg[`author-url`]}>{pkg.author.split(`<`)[0]}</a>
  (<a href="{pkg.repository}/blob/main/license">2022</a>)
</Footer>

<style>
  :global(aside.toc > nav > ol > li > a) {
    color: inherit;
  }
  :global(nav:not(.mobile) .menu) {
    gap: 1.5em;
  }
  /* On phones the fixed corner covers the top-right of the metrics table; the footer still
     links to the repo */
  @media (max-width: 600px) {
    :global(#github-corner) {
      display: none;
    }
  }
  button.find-page {
    --icon-size: 1.25em;
    display: inline-grid;
    place-items: center;
    width: 1.8em;
    height: 1.8em;
    font: inherit;
    padding: 0;
    border-radius: 50%;
    background: transparent;
  }
  button.find-page:hover {
    background: light-dark(rgba(0, 0, 100, 0.1), rgba(200, 200, 255, 0.1));
  }
  :global(::highlight(find-match)) {
    color: light-dark(#161000, #fff7c2);
    background: light-dark(#ffe066, #8a6500);
  }
  :root[data-theme='light'] img {
    filter: brightness(0.2);
  }
</style>
