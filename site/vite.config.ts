// Type declarations for Node.js modules
/// <reference types="node" />

import yaml_plugin from '@rollup/plugin-yaml'
import { sveltekit } from '@sveltejs/kit/vite'
import yaml from 'js-yaml'
import type { JSONSchema4 } from 'json-schema'
import { compile as json_to_ts } from 'json-schema-to-typescript'
import { execFile } from 'node:child_process'
import fs from 'node:fs'
import path from 'node:path'
import type { PluginOption } from 'vite'
import { defineConfig } from 'vite'

// custom Vite plugin that watches for changes to *-schema.yml files and
// automatically converts them to *-schema.d.ts files
function yaml_schema_to_typescript_plugin(): PluginOption {
  // convert *-schema.yml to *-schema.d.ts
  async function yaml_schema_to_ts(file: string): Promise<boolean> {
    try {
      const yaml_content = fs.readFileSync(file, `utf-8`)
      const file_dir = path.dirname(file)

      // Replace relative file paths in $refs with absolute file URIs
      const parsed_yaml = yaml.load(yaml_content) as JSONSchema4
      const base_name = path.basename(file, `.yml`)

      // Convert model schema to TypeScript
      const ts_name = {
        'model-schema': `ModelMetadata`,
        'dataset-schema': `DatasetRecord`,
        'label-schema': `Label`,
      }[base_name]

      const model_metadata_ts = await json_to_ts(
        parsed_yaml,
        ts_name ?? `missing ${base_name} schema`,
        {
          // should match lineWidth in deno.jsonc
          style: { semi: false, singleQuote: true, printWidth: 90 },
          bannerComment:
            `// This file is auto-generated from ${base_name}.yml. Do not edit directly.`,
          cwd: file_dir, // important for $ref resolution between *-schema.yml files
          unreachableDefinitions: true,
        },
      )

      // Write the TypeScript interface file for model schema
      const dts_file = path.resolve(`./src/lib/${base_name}.d.ts`)
      fs.writeFileSync(dts_file, model_metadata_ts)

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
      const schema_files = fs
        .readdirSync(`../tests`)
        .filter((file_name) => file_name.endsWith(`-schema.yml`))
        .map((file_name) => {
          const schema_file = path.resolve(`../tests/${file_name}`)
          yaml_schema_to_ts(schema_file)
          server.watcher.add(schema_file) // Watch each schema file for changes
          return schema_file
        })

      server.watcher.on(`change`, (changed_file) => {
        if (schema_files.includes(changed_file)) yaml_schema_to_ts(changed_file)
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
    server: {
      deps: {
        // Force Vitest to inline packages with directory imports to handle them properly
        inline: [`matterviz`, `svelte-multiselect`, `@threlte/core`, `@threlte/extras`],
      },
    },
  },

  resolve: {
    conditions: mode === `test` ? [`browser`] : undefined,
    alias: mode === `test`
      ? [
        // Mock wasm-dependent modules to avoid loading issues in jsdom
        {
          find: /^@spglib\/moyo-wasm.*/,
          replacement: new URL(`./tests/mocks/moyo-wasm.ts`, import.meta.url).pathname,
        },
      ]
      : [],
  },
}))
