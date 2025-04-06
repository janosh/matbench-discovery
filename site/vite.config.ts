// Type declarations for Node.js modules
/// <reference types="node" />

import yaml_plugin from '@rollup/plugin-yaml'
import { sveltekit } from '@sveltejs/kit/vite'
import fs from 'fs/promises'
import yaml from 'js-yaml'
import type { JSONSchema4 } from 'json-schema'
import { compile as json_to_ts } from 'json-schema-to-typescript'
import { execFile } from 'node:child_process'
import path from 'path'
import type { PluginOption } from 'vite'
import { defineConfig } from 'vite'

// custom Vite plugin that watches for changes to *-schema.yml files and
// automatically converts them to *-schema.d.ts files
function yaml_schema_to_typescript_plugin(): PluginOption {
  // convert *-schema.yml to *-schema.d.ts
  async function yaml_schema_to_ts(file: string): Promise<boolean> {
    // Read the package.json to get prettier config
    const pkg_content = await fs.readFile(path.resolve(`./package.json`), `utf-8`)
    const pkg = JSON.parse(pkg_content)
    // return if file is not a schema file

    // Read and parse the model schema YAML file
    try {
      const yaml_content = await fs.readFile(file, `utf-8`)
      const file_dir = path.dirname(file)

      // replace any relative file paths to other schema files with absolute paths
      // using the directory of the current file as the base
      const abs_path_yaml_content = yaml_content.replace(
        /\$ref: \.(.*\.yml)/g,
        (match, p1) => `\$ref: ${path.resolve(file_dir, `.${p1}`)}`,
      )
      const parsed_yaml = yaml.load(abs_path_yaml_content) as JSONSchema4
      const base_name = path.basename(file, `.yml`)

      // Convert model schema to TypeScript
      const ts_name = {
        'model-schema': `ModelMetadata`,
        'dataset-schema': `DatasetRecord`,
      }[base_name]

      const bannerComment = `// This file is auto-generated from ${base_name}.yml. Do not edit directly.`
      const model_metadata_ts = await json_to_ts(
        parsed_yaml,
        ts_name ?? `missing ${base_name} schema`,
        { style: pkg.prettier, bannerComment },
      )

      // Write the TypeScript interface file for model schema
      const dts_file = path.resolve(`./src/lib/${base_name}.d.ts`)
      await fs.writeFile(dts_file, model_metadata_ts)

      // Format model schema file
      const eslint_cmd = path.resolve(`./node_modules/.bin/eslint`)
      execFile(eslint_cmd, [`--fix`, `--config`, `eslint.config.js`, dts_file])
      return true
    } catch (error) {
      console.error(`Error processing schema file ${file}:`, error)
      return false
    }
  }

  return {
    name: `yaml-schema-to-typescript`,
    configureServer(server) {
      const model_yaml_path = path.resolve(`../tests/model-schema.yml`)
      const dataset_yaml_path = path.resolve(`../tests/dataset-schema.yml`)

      // initial TS update on server start
      yaml_schema_to_ts(model_yaml_path)
      yaml_schema_to_ts(dataset_yaml_path)

      // Watch both schema files
      server.watcher.add(model_yaml_path)
      server.watcher.add(dataset_yaml_path)

      server.watcher.on(`change`, (file) => {
        if (file.endsWith(`-schema.yml`)) yaml_schema_to_ts(file)
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
