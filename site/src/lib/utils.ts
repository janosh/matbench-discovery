export function si_fmt(num: number, fmt = `.1f`): string {
  const suffixes = [``, `k`, `M`, `B`, `T`]
  let order = Math.floor(Math.log10(Math.abs(num)) / 3)
  order = Math.min(order, suffixes.length - 1)

  const scaled = num / Math.pow(10, order * 3)
  return `${scaled.toFixed(fmt === `.0f` ? 0 : 1)}${suffixes[order]}`
}
