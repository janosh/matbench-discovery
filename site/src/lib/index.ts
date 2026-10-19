// Stringify values for display, rendering nullish/empty values as 'n/a'
export function arr_to_str(value: unknown): string {
  if (value === null || value === undefined || value === ``) return `n/a`
  if (Array.isArray(value)) return value.join(`, `)
  if (typeof value === `string`) return value
  return JSON.stringify(value)
}

// Fixed locale and UTC so prerendered HTML doesn't depend on the build machine and
// date-only strings (parsed as UTC midnight) don't shift a day west of Greenwich.
const date_formatter = new Intl.DateTimeFormat(`en-US`, {
  year: `numeric`,
  month: `short`,
  day: `numeric`,
  timeZone: `UTC`,
})
export const format_date = (date: string | number): string =>
  date_formatter.format(new Date(date))

// Compare models by benchmark inclusion date, newest first.
export const by_benchmark_added_desc = (
  model_1: { dates: { benchmark_added: string | null } },
  model_2: { dates: { benchmark_added: string | null } },
): number =>
  new Date(model_2.dates.benchmark_added ?? 0).getTime() -
  new Date(model_1.dates.benchmark_added ?? 0).getTime()
