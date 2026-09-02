<script lang="ts">
  import { goto } from '$app/navigation'
  import { page } from '$app/state'
  import { MODELS } from '$lib/models.svelte'
  import { bind_comparison_url } from '$lib/model-comparison.svelte'
  import ModelComparison from '$lib/model/ModelComparison.svelte'
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
  import type { FooterLink } from 'svelte-widgets'
  import { Changelog, Email, GitHub, RSS, Search } from 'svelte-widgets/icons'
  import MODELING_TASKS from '$pkg/modeling-tasks.yml'
  import pkg from '$site/package.json'
  import { tick, type Snippet } from 'svelte'
  import { heading_anchors } from 'svelte-widgets/heading-anchors'
  // oxlint-disable-next-line no-unassigned-import
  import '../app.css'

  let { children }: { children?: Snippet } = $props()
  let toc_desktop = $state(true)
  let find_open = $state(false)
  let find_bar = $state<ReturnType<typeof FindBar>>()
  let main_element = $state<HTMLElement>()
  bind_comparison_url()

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
  const task_labels = Object.fromEntries(
    Object.entries(MODELING_TASKS).map(([key, task]) => [
      `/tasks/${key.replaceAll(`_`, `-`)}`,
      task.label,
    ]),
  )
  // Static child pages render as dropdowns under their top-level parent route.
  const child_routes = Object.keys(import.meta.glob(`./*/*/**/+page.{svelte,md}`))
    .filter((filename) => !filename.includes(`[`))
    .map((filename) => `/${filename.split(`/`).slice(1, -1).join(`/`)}`)
    .filter((route) => !(route.startsWith(`/tasks/`) && route.includes(`/tmi`)))
  const routes = Object.keys(import.meta.glob(`./*/+page.{svelte,md}`))
    .map((filename) => `/${filename.split(`/`)[1]}`)
    .map((route) => {
      const sub_routes = child_routes.filter((child) => child.startsWith(`${route}/`))
      // include the parent route itself so Nav keeps it a clickable link (its own
      // +page) above the dropdown; Nav filters the duplicate out of the submenu
      return sub_routes.length ? { href: route, children: [route, ...sub_routes] } : route
    })
  const route_href = (route: string | { href: string }) =>
    typeof route === `string` ? route : route.href
  const nav_order: Record<string, number> = { [`/tasks`]: 0, [`/models`]: 1 }
  const ordered_routes = routes
    .filter((route) => route_href(route) !== `/changelog`)
    .toSorted(
      (route_a, route_b) =>
        (nav_order[route_href(route_a)] ?? 2) - (nav_order[route_href(route_b)] ?? 2),
    )

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
    '/data': `Details about provenance, chemistry and energies in the benchmark's train and test set.`,
    '/data/tmi': `Too much information on the benchmark's data.`,
    '/api': `API docs for the Matbench Discovery PyPI package.`,
    '/contribute': `Steps for contributing a new model to the benchmark.`,
    '/models': `Details on each model sortable by metrics.`,
    '/tasks': `Overview of all benchmark tasks for ML force fields.`,
    '/tasks/discovery': `Metrics and analysis of crystal stability prediction on the WBM test set.`,
    '/tasks/discovery/tmi': `Detailed diagnostics for the crystal discovery task.`,
    '/tasks/diatomics': `Metrics and analysis of predicting diatomic energies.`,
    '/tasks/phonons': `Metrics and analysis of predicting phonon modes and frequencies.`,
    '/tasks/geo-opt': `Metrics and analysis of predicting ground state geometries.`,
    '/tasks/md': `Metrics and analysis of molecular dynamics observables vs ab-initio reference trajectories.`,
  }
  let description = $derived(descriptions[url] ?? base_description)
  let title = $derived(url === `/` ? `` : `${url} • `)

  const actions = Object.keys(import.meta.glob(`./**/+page.{svelte,md}`))
    .filter((filename) => !filename.includes(`[`))
    .map((filename) => {
      const parts = filename.split(`/`).filter((part) => !part.startsWith(`(`)) // Remove hidden route segments
      const route = `/${parts.slice(1, -1).join(`/`)}`

      return { label: route, action: () => goto(route) }
    })
    .concat(
      MODELS.map((model) => ({
        label: `/models/${model.model_key}`,
        action: () => goto(`/models/${model.model_key}`),
      })),
    )
</script>

<CommandMenu
  {actions}
  placeholder="Go to..."
  dialog_props={{ style: `top: 15vh; bottom: auto; overflow: visible` }}
/>
<CopyButton global />

<svelte:head>
  <title>{title}Matbench Discovery</title>
  <meta name="description" content={description} />
</svelte:head>

{#if ![`/`, `/models`, `/tasks/diatomics`, `/tasks/geo-opt`].includes(url)}
  <Toc
    headingSelector={heading_selector}
    breakpoint={1350}
    minItems={3}
    hideOnIntersect="section.full-bleed .table-container, .bleed-1400"
    bind:desktop={toc_desktop}
    asideProps={{
      style: toc_desktop
        ? `max-width: 22em; position: fixed; left: calc(50vw + var(--main-max-width) / 2); top: 8em;`
        : `z-index: 1;`,
    }}
    navProps={{
      style: toc_desktop
        ? `font-size: 7pt;`
        : `font-size: 7pt; z-index: 10; padding: 1em;`,
    }}
    titleProps={{ style: `margin: 3pt` }}
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

<!-- menu_props: the svelte-widgets mobile menu hugs its content and anchors to the start
     edge, so page text shows beside the open menu; spanning the viewport fixes that.
     `max-width` clears its 90vw cap. The desktop gap lives in the style block below since an
     inline gap would also override the mobile menu's tight row spacing. -->
<Nav
  {page}
  routes={[`/`, ...ordered_routes, [pkg.paper, `Paper`]]}
  style="margin-block: 1em 0"
  menu_props={{ style: `inset-inline: 0.5rem; max-width: none` }}
  labels={{
    '/': `Home`,
    '/api': `API`,
    '/data/sets': `Datasets`,
    '/data/tmi': `TMI`,
    '/tasks/discovery/tmi': `TMI`,
    ...task_labels,
  }}
  --nav-item-padding="0 3pt"
  --nav-dropdown-link-padding="2pt 4pt"
  --nav-link-active-color="var(--link-color)"
  --nav-mobile-z-index="50"
  --nav-toggle-btn-z-index="50"
>
  {#if find_enabled}
    <button
      aria-label="Find in page"
      class="find-page"
      onclick={open_find}
      title="Find in page"
      type="button"
    >
      <Icon icon={Search} />
    </button>
  {/if}
  <ThemeToggle />
</Nav>

<main
  bind:this={main_element}
  class:bleed-1400={url === `/tasks/diatomics`}
  {@attach heading_anchors()}
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

<ModelComparison />

<Footer links={footer_links} style="--footer-bg: var(--nav-bg)">
  <img src="/favicon.svg" alt="Logo" width="30px" style="vertical-align: middle" />
  &ensp;{pkg.title} &ensp; | &ensp; ©
  <a href={pkg[`author-url`]}>{pkg.author.split(`<`)[0]}</a>
  (<a href="{pkg.repository}/blob/-/license">2022</a>)
</Footer>

<style>
  :global(aside.toc > nav > ol > li > a) {
    color: inherit;
  }
  :global(nav:not(.mobile) .menu) {
    gap: 1.5em;
  }
  /* Shim reproducing the svelte-widgets fix until it ships (see its Nav.svelte). Its mobile
     rule pads `.menu > span` — the element that paints a plain link's pill — but `.dropdown`,
     whose *child* paints the pill, so plain-link rows rendered wider and taller pills than
     dropdown rows. Move the padding onto the pill in both cases. `!important` because the
     widget's rules carry its scoping class and out-specify anything writable from here. */
  :global(nav.mobile .dropdown) {
    padding: 0 !important;
  }
  :global(nav.mobile .menu > span),
  :global(nav.mobile .dropdown > div:first-child) {
    padding: 0 6pt !important;
  }
  :global(nav.mobile .menu > span > a) {
    flex: 1;
    padding: var(--nav-item-padding, 1pt 4pt) !important;
  }
  :global(nav.mobile .dropdown > div:last-child a) {
    /* 8pt of its own indent plus the 8pt the `.dropdown` padding no longer contributes */
    margin-inline-start: 16pt !important;
    padding-block: 1pt !important;
  }
  /* Compact rows: the widget's padding is gone above, so this min-height sets row height */
  :global(nav.mobile .menu > span > a),
  :global(nav.mobile .dropdown > div:first-child > :is(a, span)),
  :global(nav.mobile .dropdown > div:last-child a) {
    display: inline-flex;
    align-items: center;
    min-height: 1.5rem;
    box-sizing: border-box;
  }
  :global(nav.mobile .dropdown > div:first-child > button) {
    min-height: 1.5rem;
    /* Finger-sized target grown leading-side only so the caret still hugs the row's
       trailing edge; the row padding is the only gap wanted */
    min-width: 2.5rem;
    justify-content: flex-end;
    padding-inline-end: 0 !important;
  }
  /* Same shim: the widget offsets the open burger's two strokes by a hardcoded 0.4rem, but
     under `space-around` adjacent bar centres are height/3 apart, so the X read lopsided. */
  :global(nav .burger[aria-expanded='true'] span:first-child) {
    transform: translateY(calc(1.4rem / 3)) rotate(45deg) !important;
  }
  :global(nav .burger[aria-expanded='true'] span:nth-child(3)) {
    transform: translateY(calc(1.4rem / -3)) rotate(-45deg) !important;
  }
  /* On phones the fixed corner covers the top-right of the metrics table; the footer still
     links to the repo */
  @media (max-width: 600px) {
    :global(#github-corner) {
      display: none;
    }
  }
  button.find-page {
    display: inline-grid;
    place-items: center;
    width: 1.8em;
    height: 1.8em;
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
