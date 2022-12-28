import yaml from '@rollup/plugin-yaml'
import { sveltekit } from '@sveltejs/kit/vite'
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
