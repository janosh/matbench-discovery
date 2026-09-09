import { pareto_staircase, sota_frontier_indices, sota_step_line } from '$lib/sota'
import { describe, expect, it } from 'vitest'

const day = (n: number) => n * 86_400_000

describe(`sota_frontier_indices`, () => {
  it.each([
    {
      desc: `strictly improving series keeps all points`,
      points: [
        { date: day(1), value: 0.5 },
        { date: day(2), value: 0.7 },
      ],
      better: `higher` as const,
      expected: [0, 1],
    },
    {
      desc: `regressions are skipped`,
      points: [
        { date: day(1), value: 0.5 },
        { date: day(2), value: 0.4 },
        { date: day(3), value: 0.6 },
      ],
      better: `higher` as const,
      expected: [0, 2],
    },
    {
      desc: `equal values don't count as new records`,
      points: [
        { date: day(1), value: 0.5 },
        { date: day(2), value: 0.5 },
      ],
      better: `higher` as const,
      expected: [0],
    },
    {
      desc: `lower-better metrics flip the comparison`,
      points: [
        { date: day(1), value: 1.2 },
        { date: day(2), value: 1.5 },
        { date: day(3), value: 0.8 },
      ],
      better: `lower` as const,
      expected: [0, 2],
    },
    {
      desc: `unsorted input is processed in date order`,
      points: [
        { date: day(3), value: 0.9 },
        { date: day(1), value: 0.5 },
        { date: day(2), value: 0.7 },
      ],
      better: `higher` as const,
      expected: [1, 2, 0],
    },
    {
      desc: `same-date tie: only the better model defines the frontier`,
      points: [
        { date: day(1), value: 0.5 },
        { date: day(2), value: 0.9 },
        { date: day(2), value: 0.7 },
      ],
      better: `higher` as const,
      expected: [0, 1],
    },
    {
      desc: `same-date tie with the worse point first in input order`,
      points: [
        { date: day(1), value: 0.5 },
        { date: day(2), value: 0.7 },
        { date: day(2), value: 0.9 },
      ],
      better: `higher` as const,
      expected: [0, 2],
    },
    {
      desc: `non-finite values never define the frontier`,
      points: [
        { date: day(1), value: NaN },
        { date: day(2), value: 0.5 },
        { date: day(3), value: Infinity },
      ],
      better: `higher` as const,
      expected: [1],
    },
    { desc: `empty input`, points: [], better: `higher` as const, expected: [] },
  ])(`$desc`, ({ points, better, expected }) => {
    expect(sota_frontier_indices(points, better)).toEqual(expected)
  })
})

describe(`sota_step_line`, () => {
  it.each([
    {
      desc: `horizontal steps with vertical jumps at each record`,
      records: [
        { date: day(1), value: 0.5 },
        { date: day(3), value: 0.8 },
      ],
      end: day(5),
      x: [day(1), day(3), day(3), day(5)],
      y: [0.5, 0.5, 0.8, 0.8],
    },
    { desc: `empty records give empty line`, records: [], end: day(9), x: [], y: [] },
    {
      desc: `single record extends to end date`,
      records: [{ date: day(2), value: 0.6 }],
      end: day(4),
      x: [day(2), day(4)],
      y: [0.6, 0.6],
    },
    {
      desc: `end date equal to last record adds no extension`,
      records: [{ date: day(2), value: 0.6 }],
      end: day(2),
      x: [day(2)],
      y: [0.6],
    },
  ])(`$desc`, ({ records, end, x, y }) => {
    expect(sota_step_line(records, end)).toEqual({ x, y })
  })
})

describe(`pareto_staircase`, () => {
  // cost-vs-fidelity shape: x lower=better (e.g. wall time), y higher=better (score)
  const points = [
    { x: 1, y: 0.5 }, // frontier: cheapest
    { x: 2, y: 0.8 }, // frontier: better but costlier
    { x: 3, y: 0.7 }, // dominated by (2, 0.8)
    { x: 4, y: 0.9 }, // frontier: best but most expensive
  ]

  it.each([
    { input: points },
    { input: [...points, points[1], { x: 2, y: 0.6 }, { x: 3, y: 0.8 }].toReversed() },
  ])(
    `keeps non-dominated points, handles ties, and inserts staircase corners`,
    ({ input }) => {
      const line = pareto_staircase(input, `lower`, `higher`)
      // corner at (next.x, current.y) between consecutive frontier points
      expect(line?.x).toEqual([1, 2, 2, 4, 4])
      expect(line?.y).toEqual([0.5, 0.5, 0.8, 0.8, 0.9])
    },
  )

  it.each([NaN, Infinity, -Infinity])(`rejects non-finite coordinates (%s)`, (value) => {
    for (const point of [
      { x: value, y: 1 },
      { x: 1, y: value },
    ]) {
      expect(() => pareto_staircase([points[0], point], `lower`, `higher`)).toThrow(
        `Pareto coordinates must be finite`,
      )
    }
  })

  it(`flips domination with axis directions`, () => {
    // (4, 0.9) has the best x and (1, 0.5) the best y; both remain on the frontier.
    const line = pareto_staircase(points, `higher`, `lower`)
    expect(line?.x[0]).toBe(4) // sorted best-x first under higher-is-better x
    expect(line?.y.at(-1)).toBe(0.5) // ends at the best-y point
  })

  it.each([1, 200_000])(
    `draws an L through an all-dominating point (%s repeated cohorts)`,
    (repetitions) => {
      const pts = [
        { x: 1, y: 0.9 },
        { x: 2, y: 0.5 },
        { x: 3, y: 0.1 },
      ]
      const line = pareto_staircase(
        Array.from({ length: repetitions }, () => pts).flat(),
        `lower`,
        `higher`,
      )
      // vertical lead-in from the worst-y extent, horizontal tail-out to worst-x extent
      expect(line?.x).toEqual([1, 1, 3])
      expect(line?.y).toEqual([0.1, 0.9, 0.9])
    },
  )

  it.each([
    [`empty input`, []],
    [`single point (no region to frame)`, [{ x: 1, y: 1 }]],
  ])(`returns null for %s`, (_name, pts) => {
    expect(pareto_staircase(pts, `lower`, `higher`)).toBeNull()
  })
})
