// Running-best ("SOTA") frontier helpers for the model progress timeline.

export interface SotaPoint {
  date: number // unix epoch ms
  value: number
}

// Indices (into the input array) of points that set a new record at their date.
// Points are processed in date order; same-date ties are processed best-first so
// only the best model of a given day can define the frontier. Only strict
// improvements count as new records.
export const sota_frontier_indices = (
  points: readonly SotaPoint[],
  better: `higher` | `lower` = `higher`,
): number[] => {
  const sign = better === `higher` ? 1 : -1
  const order = points
    .map((_, idx) => idx)
    .toSorted(
      (idx_a, idx_b) =>
        points[idx_a].date - points[idx_b].date ||
        sign * (points[idx_b].value - points[idx_a].value),
    )
  const records: number[] = []
  let best = -Infinity
  for (const idx of order) {
    const score = sign * points[idx].value
    if (Number.isFinite(score) && score > best) {
      best = score
      records.push(idx)
    }
  }
  return records
}

// Step-line coordinates tracing the running best: horizontal segments at each
// record's value, jumping vertically at the date of the next record, extended to
// `end_date` (typically today) so the current SOTA reads as "still standing".
export const sota_step_line = (
  records: readonly SotaPoint[],
  end_date: number,
): { x: number[]; y: number[] } => {
  const x: number[] = []
  const y: number[] = []
  for (const [idx, record] of records.entries()) {
    if (idx > 0) {
      // horizontal run-up to this record's date at the previous value
      x.push(record.date)
      y.push(records[idx - 1].value)
    }
    x.push(record.date)
    y.push(record.value)
  }
  if (records.length > 0 && end_date > records[records.length - 1].date) {
    x.push(end_date)
    y.push(records[records.length - 1].value)
  }
  return { x, y }
}
