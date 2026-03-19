// Type declarations for Node.js modules
/// <reference types="node" />

import yaml_plugin from '@rollup/plugin-yaml'
import { sveltekit } from '@sveltejs/kit/vite'
import yaml from 'js-yaml'
import type { JSONSchema4 } from 'json-schema'
import { compile as json_to_ts } from 'json-schema-to-typescript'
import { execFileSync } from 'node:child_process'
import fs from 'node:fs'
import path from 'node:path'
import type { Plugin } from 'vite'
import { defineConfig } from 'vite-plus'

// custom Vite plugin that watches for changes to *-schema.yml files and
// automatically converts them to *-schema.d.ts files
function yaml_schema_to_typescript_plugin(): Plugin {
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
          // should match fmt config in vite.config.ts
          style: { semi: false, singleQuote: true, printWidth: 90 },
          bannerComment: `// This file is auto-generated from ${base_name}.yml. Do not edit directly.`,
          cwd: file_dir, // important for $ref resolution between *-schema.yml files
          unreachableDefinitions: true,
        },
      )

      // Write the TypeScript interface file for model schema
      const dts_file = path.resolve(`./src/lib/${base_name}.d.ts`)
      fs.writeFileSync(dts_file, model_metadata_ts)

      // Format generated schema file
      const vp_cmd = path.resolve(`./node_modules/.bin/vp`)
      execFileSync(vp_cmd, [`fmt`, `--write`, dts_file])
      return true
    } catch (error) {
      console.error(`Error processing schema file ${file}:`, error)
      return false
    }
  }

  return {
    name: `yaml-schema-to-typescript`,
    configureServer(server) {
      const schema_files = new Set(
        fs
          .readdirSync(`../tests`)
          .filter((file_name) => file_name.endsWith(`-schema.yml`))
          .map((file_name) => {
            const schema_file = path.resolve(`../tests/${file_name}`)
            void yaml_schema_to_ts(schema_file)
            server.watcher.add(schema_file) // Watch each schema file for changes
            return schema_file
          }),
      )

      server.watcher.on(`change`, (changed_file) => {
        if (schema_files.has(changed_file)) void yaml_schema_to_ts(changed_file)
      })
    },
  }
}

export default defineConfig({
  fmt: {
    semi: false,
    singleQuote: true,
    printWidth: 90,
    ignorePatterns: [
      `src/figs/**/*.svelte`,
      `src/figs/**/*.json`,
      `src/routes/**/*.json`,
    ],
  },
  lint: {
    plugins: [`oxc`, `typescript`, `unicorn`, `import`, `jest`],
    options: {
      typeAware: true,
      typeCheck: true,
    },
    categories: {
      correctness: `error`,
      suspicious: `error`,
      pedantic: `error`,
      perf: `error`,
    },
    ignorePatterns: [
      `build/`,
      `.svelte-kit/`,
      `src/figs/`,
      `dist/`,
      `src/lib/*.d.ts`,
      `scripts/`,
    ],
    rules: {
      '@typescript-eslint/no-explicit-any': `error`,
      '@typescript-eslint/no-non-null-asserted-optional-chain': `error`,
      '@typescript-eslint/no-non-null-assertion': `error`,
      'no-unused-vars': `off`, // superseded by type-aware version below
      '@typescript-eslint/no-unused-vars': [
        `error`,
        { argsIgnorePattern: `^_`, varsIgnorePattern: `^_` },
      ],
      'no-eval': `error`,
      eqeqeq: `error`,
      'no-var': `error`,
      'no-throw-literal': `error`,
      'no-useless-rename': `error`,
      'no-self-compare': `error`,
      'no-template-curly-in-string': `error`,
      'no-constructor-return': `error`,
      'no-console': [`error`, { allow: [`warn`, `error`] }],
      'default-param-last': `error`,
      'guard-for-in': `error`,
      'require-await': `error`,
      'eslint-plugin-unicorn/no-useless-spread': `error`,
      'eslint-plugin-unicorn/prefer-string-replace-all': `error`,
      'eslint-plugin-unicorn/catch-error-name': `error`,
      'eslint-plugin-unicorn/prefer-set-has': `error`,
      'eslint-plugin-unicorn/prefer-dom-node-append': `error`,
      'eslint-plugin-import/no-duplicates': `error`,
      'no-inner-declarations': `error`,
      'eslint-plugin-unicorn/prefer-global-this': `error`,
      'eslint-plugin-unicorn/no-lonely-if': `error`,
      'eslint-plugin-unicorn/no-negated-condition': `error`,
      'eslint-plugin-unicorn/no-typeof-undefined': `error`,
      'eslint-plugin-unicorn/prefer-optional-catch-binding': `error`,
      'eslint-plugin-unicorn/no-length-as-slice-end': `error`,
      'eslint-plugin-unicorn/prefer-node-protocol': `error`,
      'eslint-plugin-unicorn/prefer-regexp-test': `error`,
      'eslint-plugin-unicorn/throw-new-error': `error`,
      'eslint-plugin-unicorn/prefer-includes': `error`,
      'eslint-plugin-unicorn/prefer-type-error': `error`,
      'eslint-plugin-unicorn/prefer-date-now': `error`,
      'eslint-plugin-unicorn/require-number-to-fixed-digits-argument': `error`,
      'eslint-plugin-unicorn/no-useless-promise-resolve-reject': `error`,
      'no-self-assign': `off`, // Svelte reactive `x = x` assignments
      'no-await-in-loop': `off`, // test code uses sequential await in loops
      'no-shadow': `off`, // closures intentionally shadow outer names
      'prefer-const': `off`, // `let` needed for $state/$derived/$bindable
      'no-negated-condition': `off`, // `!== undefined` ternaries are idiomatic
      'no-unassigned-vars': `off`, // Svelte bind: variables assigned by framework
      'no-map-spread': `off`, // functional `.map(x => ({...x, prop}))` is standard
      '@typescript-eslint/no-unnecessary-condition': `off`, // reactive narrowing false positives
      '@typescript-eslint/consistent-type-imports': `off`, // template component import false positives
      'eslint-plugin-unicorn/consistent-function-scoping': `off`, // test helpers + Svelte reactive closures
      'eslint-plugin-unicorn/no-new-array': `off`, // `new Array(n).fill()` is standard
      'eslint-plugin-unicorn/prefer-array-find': `off`, // false positives with .filter() for all matches
      'eslint-plugin-import/no-self-import': `off`, // self-mounting components
      'eslint-plugin-import/no-unassigned-import': `off`, // CSS imports are side-effect-only
      'eslint-plugin-import/no-named-as-default-member': `off`, // yaml.load() is the documented API
      // DOM/any propagation — oxlint lacks DOM type stubs
      '@typescript-eslint/no-unsafe-argument': `off`,
      '@typescript-eslint/no-unsafe-assignment': `off`,
      '@typescript-eslint/no-unsafe-call': `off`,
      '@typescript-eslint/no-unsafe-member-access': `off`,
      '@typescript-eslint/no-unsafe-return': `off`,
      // Pedantic rules too noisy for this codebase
      'no-inline-comments': `off`,
      'no-confusing-void-expression': `off`,
      'no-promise-executor-return': `off`,
      'strict-boolean-expressions': `off`, // truthiness checks are idiomatic
      'max-lines-per-function': `off`,
      'max-lines': `off`,
      'max-depth': `off`,
      'max-classes-per-file': `off`,
      'sort-vars': `off`,
      'eslint-plugin-jest/no-conditional-in-test': `off`, // parameterized tests use conditionals
      'eslint-plugin-jest/no-conditional-expect': `off`, // vitest conditional asserts in parameterized tests
      'eslint-plugin-jest/no-disabled-tests': `off`, // intentionally skipped tests pending upstream fixes
      'eslint-plugin-jest/valid-expect': `off`, // vitest supports message argument (jest rule false positive)
      'no-warning-comments': `off`, // TODO comments are legitimate development markers
      'no-else-return': `off`, // early return style is not always clearer
      '@typescript-eslint/no-unsafe-type-assertion': `off`, // `as Type` casts are standard TS practice
      '@typescript-eslint/only-throw-error': `off`, // SvelteKit redirect() throws non-Error objects
      '@typescript-eslint/no-misused-promises': `off`, // event handler async patterns
      '@typescript-eslint/no-deprecated': `off`, // false positives on standard DOM APIs (createElement)
      '@typescript-eslint/await-thenable': `off`, // sync-returning functions may be awaited for test clarity
      'eslint-plugin-unicorn/no-array-callback-reference': `off`,
      'eslint-plugin-unicorn/no-useless-undefined': `off`,
      'eslint-plugin-unicorn/no-object-as-default-parameter': `off`,
      'eslint-plugin-import/max-dependencies': `off`,
    },
  },
  staged: {
    '*.{js,ts,svelte,html,css,md,json,yaml}': `vp check --fix`,
    '*.{ts,svelte}': `sh -c 'npx svelte-kit sync && npx svelte-check-rs --threshold error'`,
    '*.test.ts': `sh -c '! grep -E "(test|describe)\\.only\\(" "$@"' --`,
    '*': `codespell --ignore-words-list falsy --check-filenames`,
  },
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
    environment: `happy-dom`, // Faster than jsdom
    css: true,
    pool: `threads`, // Parallel test execution
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
    conditions: process.env.VITEST ? [`browser`] : undefined,
    alias: process.env.VITEST
      ? [
          // Mock wasm-dependent modules to avoid loading issues in jsdom
          {
            find: /^@spglib\/moyo-wasm.*/,
            replacement: new URL(`./tests/mocks/moyo-wasm.ts`, import.meta.url).pathname,
          },
        ]
      : [],
  },
})
