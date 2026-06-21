import type { Config } from '@sveltejs/kit'
import adapter from '@sveltejs/adapter-static'
import { mdsvex } from 'mdsvex'
import pkg from './package.json' with { type: 'json' }
import katex from 'rehype-katex-svelte'
import math from 'remark-math' // Remark-math@3.0.0 pinned due to mdsvex, see https://github.com/kwshi/rehype-katex-svelte#usage
import { heading_ids } from 'svelte-multiselect/heading-anchors' // Adds IDs to headings at build time

export default {
  extensions: [`.svelte`, `.svx`, `.md`, `.html`],

  preprocess: [
    // Replace readme links to docs with site-internal links
    // (which don't require browser navigation)
    {
      markup: ({ content }) => ({ code: content.replaceAll(pkg.homepage, ``) }),
    },
    mdsvex({
      // cast bridges katex's unified version to the one mdsvex bundles
      rehypePlugins: [katex] as NonNullable<
        Parameters<typeof mdsvex>[0]
      >[`rehypePlugins`],
      remarkPlugins: [math],
      extensions: [`.svx`, `.md`],
    }),
    heading_ids(), // Runs after mdsvex converts markdown to HTML
    {
      markup: (file) => {
        const filename = file.filename?.replaceAll(`\\`, `/`) ?? ``
        const route = filename.split(`site/src/routes/`)[1] ?? ``

        // Only manuscript figure markdown uses @label/@fig citation-style rewrites.
        if (!route.includes(`discovery-metric-figs`)) return { code: file.content }

        const fig_index: string[] = []

        // Replace figure labels with 'Fig. {n}' and add to fig_index
        let code = file.content.replaceAll(
          /@label:(?<id>(?:fig|tab):[^\s]+)/g,
          (_match, id) => {
            if (!fig_index.includes(id)) fig_index.push(id)
            const idx = (route.startsWith(`si`) ? `S` : ``) + fig_index.length
            const link_icon = `<a aria-hidden="true" tabindex="-1" href="#${id}"><svg width="16" height="16" viewBox="0 0 16 16"><use xlink:href="#octicon-link"></use></svg></a>`
            return `<strong id='${id}'>${link_icon}Fig. ${idx}</strong>`
          },
        )

        // Replace figure references @fig:label with 'fig. {n}' and add to fig_index
        code = code.replaceAll(
          /@(?<id>(?<fig>fig):(?:[a-z0-9]+-?)+)/gi, // Match case-insensitive but replace case-sensitive
          // @(f|F)ig becomes '(f|F)ig. {n}'
          (_full_str, id, fig_or_Fig) => {
            const id_lower = id.toLowerCase()
            let idx: number | string = fig_index.indexOf(id_lower) + 1
            if (idx === 0) {
              console.warn(
                `Figure id='${id}' not found, expected one of ${fig_index.join(`, `)}`,
              )
              idx = `not found`
            }
            return `<a href="#${id_lower}">${fig_or_Fig}. ${idx}</a>`
          },
        )

        // Preprocess markdown citations @auth_1st-word-title_yyyy into citation links
        // Links to bibliography items, href must match id format in References.svelte
        code = code.replaceAll(
          /\[?@(?<id>(?<author>.+?)_.+?_(?<year>\d{4}));?\]?/g, // Ends with ;?\]? to match single and multiple citations
          (_match, id, author, year) =>
            `[<a class="ref" href="#${id}">${author} ${year}</a>]`,
        )

        return { code }
      },
    },
  ],

  kit: {
    adapter: adapter(),

    alias: {
      $site: `.`,
      $root: `..`,
      $data: `../data`,
      $pkg: `../matbench_discovery`,
      $figs: `src/figs`,
      $routes: `src/routes`,
    },
  },
} satisfies Config
