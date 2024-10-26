import adapter from '@sveltejs/adapter-static'
import { s } from 'hastscript'
import { mdsvex } from 'mdsvex'
import link_headings from 'rehype-autolink-headings'
import katex from 'rehype-katex-svelte'
import heading_slugs from 'rehype-slug'
import math from 'remark-math'
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
      rehypePlugins: [
        katex,
        heading_slugs,
        [
          link_headings,
          {
            behavior: `append`,
            test: [`h2`, `h3`, `h4`, `h5`, `h6`], // don't auto-link <h1>
            content: s(
              `svg`,
              { width: 16, height: 16, viewBox: `0 0 16 16` },
              // symbol #octicon-link defined in app.html
              s(`use`, { 'xlink:href': `#octicon-link` }),
            ),
          },
        ],
      ],
      // remark-math@3.0.0 pinned due to mdsvex, see
      // https://github.com/kwshi/rehype-katex-svelte#usage
      remarkPlugins: [math],
      extensions: [`.svx`, `.md`],
    }),
    {
      markup: (file) => {
        const route = file.filename.split(`site/src/routes/`)[1]?.split(`/`)[0]
        if (!route) return

        if ([`paper`, `preprint`, `si`].some((key) => route.startsWith(key))) {
          let fig_index = []
          let ref_index = []

          // Replace figure labels with 'Fig. {n}' and add to fig_index
          let code = file.content.replace(
            /@label:((fig|tab):[^\s]+)/g,
            (_match, id) => {
              if (!fig_index.includes(id)) fig_index.push(id)
              const idx = (route.startsWith(`si`) ? `S` : ``) + fig_index.length
              const link_icon = `<a aria-hidden="true" tabindex="-1" href="#${id}"><svg width="16" height="16" viewBox="0 0 16 16"><use xlink:href="#octicon-link"></use></svg></a>`
              return `<strong id='${id}'>${link_icon}Fig. ${idx}</strong>`
            },
          )

          // Replace figure references @fig:label with 'fig. {n}' and add to fig_index
          code = code.replace(
            /@((fig):([a-z0-9]+-?)+)/gi, // match case-insensitive but replace case-sensitive
            // @(f|F)ig becomes '(f|F)ig. {n}'
            (_full_str, id, fig_or_Fig) => {
              const id_lower = id.toLowerCase()
              let idx = fig_index.indexOf(id_lower) + 1
              if (idx == 0) {
                console.trace(
                  `Figure id='${id}' not found, expected one of ${fig_index}`,
                )
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
              if (idx == -1) {
                ref_index.push(id)
                idx = ref_index.length - 1
              }
              return `[<a class="ref" href="#${id}">${author} ${year}</a>]`
            },
          )

          return { code }
        }
      },
    },
  ],

  kit: {
    adapter: adapter(),

    alias: {
      $site: `.`,
      $root: `..`,
      $pkg: `../matbench_discovery`,
      $figs: `src/figs`,
    },
  },
}
