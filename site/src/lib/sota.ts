// Running-best ("SOTA") frontier helpers for the model progress timeline.

// Indices (into the input array) of points that set a new record at their date.
// Points are processed in date order; same-date ties are processed best-first so
// only the best model of a given day can define the frontier. Only strict
// improvements count as new records.
export const sota_frontier_indices = (
  points: readonly { date: number; value: number }[],
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
  records: readonly { date: number; value: number }[],
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

// Non-dominated staircase for two directed axes: keeps points no other point beats
// (at least as good on both axes, strictly better on one), sorted along x, with a
// corner inserted between consecutive frontier points so the line steps along the
// dominated-region boundary instead of cutting across. Both ends extend to the data
// extents so the boundary frames the whole dominated region - crucially, a single
// model dominating the entire field still draws an L through its point instead of
// nothing. Null only when < 2 input points (no region to frame).
export const pareto_staircase = (
  points: readonly { x: number; y: number }[],
  x_better: `higher` | `lower`,
  y_better: `higher` | `lower`,
): { x: number[]; y: number[] } | null => {
  if (points.length < 2) return null
  // sign-flip so both axes minimize
  const [sx, sy] = [x_better === `higher` ? -1 : 1, y_better === `higher` ? -1 : 1]
  const frontier = points
    .filter(
      (pt) =>
        !points.some(
          (other) =>
            sx * other.x <= sx * pt.x &&
            sy * other.y <= sy * pt.y &&
            (sx * other.x < sx * pt.x || sy * other.y < sy * pt.y),
        ),
    )
    .toSorted((p1, p2) => sx * (p1.x - p2.x))
  if (frontier.length === 0) return null // all points non-finite/mutually dominated

  const [x, y] = [[] as number[], [] as number[]]
  const push = (px: number, py: number) => {
    if (px === x.at(-1) && py === y.at(-1)) return // skip zero-length segments
    x.push(px)
    y.push(py)
  }
  // data extremes in the "worse" direction of each axis
  const worst_x = (sx > 0 ? Math.max : Math.min)(...points.map((pt) => pt.x))
  const worst_y = (sy > 0 ? Math.max : Math.min)(...points.map((pt) => pt.y))

  push(frontier[0].x, worst_y) // vertical lead-in at the best-x end
  push(frontier[0].x, frontier[0].y)
  for (const [idx, pt] of frontier.slice(1).entries()) {
    push(pt.x, frontier[idx].y) // corner along the dominated-region boundary
    push(pt.x, pt.y)
  }
  push(worst_x, frontier[frontier.length - 1].y) // horizontal tail-out
  return x.length < 2 ? null : { x, y }
}
