import { sota_frontier_indices, sota_step_line } from '$lib/sota'
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
