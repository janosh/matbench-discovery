// Mock for @spglib/moyo-wasm to avoid wasm loading issues in jsdom tests

// Default export for init function
export default async function init(): Promise<void> {}

export class MoyoDataset {
  constructor() {}
}

export const get_version = (): string => `0.0.0-mock`

export const analyze_cell = (): null => null

// Mock wasm url import
export const moyo_wasm_url = ``
