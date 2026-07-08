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
  it(`builds horizontal steps with vertical jumps at each record`, () => {
    const records = [
      { date: day(1), value: 0.5 },
      { date: day(3), value: 0.8 },
    ]
    const { x, y } = sota_step_line(records, day(5))
    expect(x).toEqual([day(1), day(3), day(3), day(5)])
    expect(y).toEqual([0.5, 0.5, 0.8, 0.8])
  })

  it.each([
    { desc: `empty records give empty line`, records: [], end: day(9), x: [], y: [] },
    {
      desc: `single record extends to end date`,
      records: [{ date: day(2), value: 0.6 }],
      end: day(4),
      x: [day(2), day(4)],
      y: [0.6, 0.6],
    },
    {
      desc: `end date before last record adds no extension`,
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

  it(`keeps non-dominated points and inserts staircase corners`, () => {
    const line = pareto_staircase(points, `lower`, `higher`)
    // corner at (next.x, current.y) between consecutive frontier points
    expect(line?.x).toEqual([1, 2, 2, 4, 4])
    expect(line?.y).toEqual([0.5, 0.5, 0.8, 0.8, 0.9])
  })

  it(`flips domination with axis directions`, () => {
    // x higher=better, y lower=better: (4, 0.9) is now worst on both -> dominated by
    // any point with lower y... actually only frontier points survive domination
    const line = pareto_staircase(points, `higher`, `lower`)
    // frontier sorted by descending x: (4, 0.9) dominated by (2,0.8)? x=4 better,
    // y=0.9 worse -> not dominated. (1, 0.5) has best y but worst x -> frontier.
    expect(line?.x[0]).toBe(4) // sorted best-x first under higher-is-better x
    expect(line?.y.at(-1)).toBe(0.5) // ends at the best-y point
  })

  it(`draws an L through a single all-dominating point (e.g. fastest AND best model)`, () => {
    const pts = [
      { x: 1, y: 0.9 },
      { x: 2, y: 0.5 },
      { x: 3, y: 0.1 },
    ]
    const line = pareto_staircase(pts, `lower`, `higher`)
    // vertical lead-in from the worst-y extent, horizontal tail-out to worst-x extent
    expect(line?.x).toEqual([1, 1, 3])
    expect(line?.y).toEqual([0.1, 0.9, 0.9])
  })

  it.each([
    [`empty input`, []],
    [`single point (no region to frame)`, [{ x: 1, y: 1 }]],
  ])(`returns null for %s`, (_name, pts) => {
    expect(pareto_staircase(pts, `lower`, `higher`)).toBeNull()
  })
})
