// Stringify values for display, rendering nullish/empty values as 'n/a'
export function arr_to_str(value: unknown): string {
  if (value === null || value === undefined || value === ``) return `n/a`
  if (Array.isArray(value)) return value.join(`, `)
  return JSON.stringify(value)
}

export const format_date = (
  date: string | number,
  options?: Intl.DateTimeFormatOptions,
): string =>
  new Date(date).toLocaleDateString(undefined, {
    year: `numeric`,
    month: `short`,
    day: `numeric`,
    ...options,
  })

// Compare models by benchmark inclusion date, newest first.
export const by_benchmark_added_desc = (
  model_1: { dates: { benchmark_added: string | null } },
  model_2: { dates: { benchmark_added: string | null } },
): number =>
  new Date(model_2.dates.benchmark_added ?? 0).getTime() -
  new Date(model_1.dates.benchmark_added ?? 0).getTime()
