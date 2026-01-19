import adapter from '@sveltejs/adapter-static'
import { mdsvex } from 'mdsvex'
import katex from 'rehype-katex-svelte'
import math from 'remark-math' // remark-math@3.0.0 pinned due to mdsvex, see https://github.com/kwshi/rehype-katex-svelte#usage
import { heading_ids } from 'svelte-multiselect/heading-anchors' // adds IDs to headings at build time
import { sveltePreprocess } from 'svelte-preprocess'

const { default: pkg } = await import(`./package.json`, {
  with: { type: `json` },
})

/** @type {import('@sveltejs/kit').Config} */
export default {
  extensions: [`.svelte`, `.svx`, `.md`, `.html`],

  preprocess: [
    // replace readme links to docs with site-internal links
    // (which don't require browser navigation)
    sveltePreprocess({ replace: [[pkg.homepage, ``]] }),
    mdsvex({
      rehypePlugins: [katex],
      remarkPlugins: [math],
      extensions: [`.svx`, `.md`],
    }),
    heading_ids(), // runs after mdsvex converts markdown to HTML
    {
      markup: (file) => {
        const route = file.filename.split(`site/src/routes/`)[1]

        if (!route?.includes(`discovery-metric-figs`)) return

        const fig_index = []
        const ref_index = []

        // Replace figure labels with 'Fig. {n}' and add to fig_index
        let code = file.content.replace(/@label:((fig|tab):[^\s]+)/g, (_match, id) => {
          if (!fig_index.includes(id)) fig_index.push(id)
          const idx = (route.startsWith(`si`) ? `S` : ``) + fig_index.length
          const link_icon =
            `<a aria-hidden="true" tabindex="-1" href="#${id}"><svg width="16" height="16" viewBox="0 0 16 16"><use xlink:href="#octicon-link"></use></svg></a>`
          return `<strong id='${id}'>${link_icon}Fig. ${idx}</strong>`
        })

        // Replace figure references @fig:label with 'fig. {n}' and add to fig_index
        code = code.replace(
          /@((fig):([a-z0-9]+-?)+)/gi, // match case-insensitive but replace case-sensitive
          // @(f|F)ig becomes '(f|F)ig. {n}'
          (_full_str, id, fig_or_Fig) => {
            const id_lower = id.toLowerCase()
            let idx = fig_index.indexOf(id_lower) + 1
            if (idx === 0) {
              console.trace(`Figure id='${id}' not found, expected one of ${fig_index}`)
              idx = `not found`
            }
            return `<a href="#${id_lower}">${fig_or_Fig}. ${idx}</a>`
          },
        )

        // preprocess markdown citations @auth_1st-word-title_yyyy into superscript
        // links to bibliography items, href must match id format in References.svelte
        code = code.replace(
          /\[?@((.+?)_.+?_(\d{4}));?\]?/g, // ends with ;?\]? to match single and multiple citations
          (_match, id, author, year) => {
            let idx = ref_index.indexOf(id)
            if (idx === -1) {
              ref_index.push(id)
              idx = ref_index.length - 1
            }
            return `[<a class="ref" href="#${id}">${author} ${year}</a>]`
          },
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
}
