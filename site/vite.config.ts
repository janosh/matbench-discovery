import yaml from '@rollup/plugin-yaml'
import { sveltekit } from '@sveltejs/kit/vite'
import { exec } from 'child_process'
import { resolve } from 'path'
import type { UserConfig } from 'vite'

const vite_config: UserConfig = {
  plugins: [sveltekit(), yaml()],

  resolve: {
    alias: {
      $src: resolve(`./src`),
      $site: resolve(`.`),
      $root: resolve(`..`),
    },
  },

  server: {
    fs: { allow: [`../..`] }, // needed to import readme.md
    port: 3000,
  },

  preview: {
    port: 3000,
  },
}

export default vite_config

if (process.env.PROD) {
  // update generated API docs on production builds
  const src_url = `https://github.com/janosh/matbench-discovery/blob/main`
  const route = `src/routes/api`
  await exec(`rm -f ${route}/*.md`)
  await exec(
    `cd .. && lazydocs matbench_discovery --output-path site/${route} --no-watermark --src-base-url ${src_url}`
  )

  // remove <b> tags from generated markdown
  await exec(`sed -i 's/<b>//g' ${route}/*.md`)
  await exec(`sed -i 's/<\\/b>//g' ${route}/*.md`)
  // tweak look of badges linking to source code
  const old_src = `src="https://img.shields.io/badge/-source-cccccc?style=flat-square"`
  const new_src = `src="https://img.shields.io/badge/source-blue?style=flat" alt="source link"`
  await exec(`sed -i 's/${old_src}/${new_src}/g' ${route}/*.md`)
}
