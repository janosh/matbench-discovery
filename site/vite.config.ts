// Type declarations for Node.js modules
/// <reference types="node" />

import yaml_plugin from '@rollup/plugin-yaml'
import { sveltekit } from '@sveltejs/kit/vite'
import fs from 'fs/promises'
import yaml from 'js-yaml'
import { compile as json_to_ts } from 'json-schema-to-typescript'
import path from 'path'
import type { PluginOption } from 'vite'
import { defineConfig } from 'vite'

// custom Vite plugin that watches for changes to model-schema.yml and automatically
// regenerates the TypeScript definitions in model-schema.d.ts
function yaml_schema_to_typescript_plugin(): PluginOption {
  const yaml_schema_path = path.resolve(`../tests/model-schema.yml`)

  // convert model-schema.yml to model-schema.d.ts
  async function generate_typescript_from_schema(): Promise<boolean> {
    // Read the package.json to get prettier config
    const pkg_content = await fs.readFile(path.resolve(`./package.json`), `utf-8`)
    const pkg = JSON.parse(pkg_content)

    // Read and parse the YAML file
    const yaml_content = await fs.readFile(yaml_schema_path, `utf-8`)
    const parsed_yaml = yaml.load(yaml_content)

    // Convert schema to TypeScript
    const model_metadata_ts = await json_to_ts(parsed_yaml, `ModelMetadata`, {
      style: pkg.prettier,
      bannerComment: `// This file is auto-generated from model-schema.yml. Do not edit directly.`,
    })

    // Write the TypeScript interface file
    const dts_out_file = path.resolve(`./src/lib/model-schema.d.ts`)
    await fs.writeFile(dts_out_file, model_metadata_ts)
    return true
  }

  return {
    name: `yaml-schema-to-typescript`,
    configureServer(server) {
      // Initial model-schema.yml to model-schema.d.ts conversion when server starts
      generate_typescript_from_schema().catch((err) => {
        console.error(`Failed to generate TypeScript on startup: ${err}`)
      })
      server.watcher.add(yaml_schema_path) // watch model-schema.yml

      server.watcher.on(`change`, async (file) => {
        if (file.endsWith(`model-schema.yml`)) await generate_typescript_from_schema()
      })
    },
  }
}

export default defineConfig(({ mode }) => ({
  plugins: [
    sveltekit(),
    yaml_plugin({ extensions: [`.yml`, `.yaml`, `.cff`] }),
    yaml_schema_to_typescript_plugin(),
  ],

  server: {
    fs: { allow: [`../..`] }, // needed to import from $root
    port: 3000,
  },

  preview: {
    port: 3000,
  },

  test: {
    environment: `jsdom`,
    css: true,
    coverage: {
      reporter: [`text`, `json-summary`],
    },
    setupFiles: `tests/index.ts`,
    dir: `tests`,
    include: [`**/*.test.ts`, `**/*.test.svelte.ts`],
  },

  resolve: {
    conditions: mode === `test` ? [`browser`] : [],
  },
}))
