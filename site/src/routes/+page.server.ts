import schema from '$root/tests/model-metadata-schema.yaml'
import fs from 'fs'
import { compile } from 'json-schema-to-typescript'

// compile TypeScript types from YAML schema
const model_md_type = await compile(schema, `ModelMetadata`, {
  style: { singleQuote: true, semi: false, printWidth: 100 },
})
const dts_out_file = `src/lib/model-metadata.d.ts`
fs.writeFileSync(dts_out_file, model_md_type)
