// Mock for @spglib/moyo-wasm to avoid wasm loading issues in jsdom tests

// Default export for init function
export default async function init(): Promise<void> {}

export class MoyoDataset {
  name = `mock`
}

export const get_version = (): string => `0.0.0-mock`

export const analyze_cell = (): null => null

// space-group database helpers matterviz imports from wyckoff-db.js
export const hall_symbol_entries_from_number = (): unknown[] => []

export const wyckoff_positions = (): unknown[] => []

// Mock wasm url import
export const moyo_wasm_url = ``
