import adapter from '@sveltejs/adapter-static'
import { s } from 'hastscript'
import { mdsvex } from 'mdsvex'
import link_headings from 'rehype-autolink-headings'
import katex from 'rehype-katex-svelte'
import heading_slugs from 'rehype-slug'
import math from 'remark-math'
import preprocess from 'svelte-preprocess'
import assets from 'svelte-preprocess-import-assets'

const rehypePlugins = [
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
        s(`use`, { 'xlink:href': `#octicon-link` })
      ),
    },
  ],
]

const { default: pkg } = await import(`./package.json`, {
  assert: { type: `json` },
})

/** @type {import('@sveltejs/kit').Config} */
export default {
  extensions: [`.svelte`, `.svx`, `.md`],

  preprocess: [
    {
      markup: (file) => {
        if (file.filename.endsWith(`paper/+page.svx`)) {
          // preprocess markdown citations @auth_1st-word-title_yyyy into superscript
          // links to bibliography items, href must match id format in References.svelte
          const code = file.content.replace(
            /@((.+?)_.+?_(\d{4}))/g,
            (_full_str, bib_id, author, year) =>
              `<sup><a href="#${bib_id}">${author} ${year}</a></sup>`
          )
          return { code }
        }
      },
    },
    // replace readme links to docs with site-internal links
    // (which don't require browser navigation)
    preprocess({ replace: [[pkg.homepage, ``]] }),
    mdsvex({
      rehypePlugins,
      // remark-math@3.0.0 pinned due to mdsvex, see
      // https://github.com/kwshi/rehype-katex-svelte#usage
      remarkPlugins: [math],
      extensions: [`.svx`, `.md`],
    }),
    assets(),
  ],

  kit: {
    adapter: adapter(),

    alias: {
      $site: `.`,
      $root: `..`,
      $static: `./static`,
      $figs: `./static/figs`,
    },
  },
}
