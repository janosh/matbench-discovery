import NodeUtility from 'happy-dom/lib/nodes/node/NodeUtility.js'
import { expect, it } from 'vite-plus/test'

type HappyDomNode = Parameters<typeof NodeUtility.isFollowing>[0]

// happy-dom 20's original isFollowing verbatim (tests/index.ts replaces it)
const walk_is_following = (node_a: HappyDomNode, node_b: HappyDomNode): boolean => {
  if (node_a === node_b) return false
  let current: HappyDomNode | null = node_b
  while (current) {
    current = NodeUtility.following(current)
    if (current === node_a) return true
  }
  return false
}

it(`patched happy-dom isFollowing matches the original document-order walk`, () => {
  expect(String(NodeUtility.isFollowing)).not.toContain(`following(`) // patch installed
  let seed = 0x2f_6b_3a_1d // deterministic LCG so failures reproduce
  const rand_int = (max: number): number => {
    seed = (Math.imul(seed, 1_664_525) + 1_013_904_223) >>> 0
    return seed % max
  }
  // random trees mixing proxied (form, select) elements, text, comments and a shadow root
  const tags = [`div`, `span`, `form`, `select`, `option`, `table`, `tr`, `td`, `a`]
  const make_tree = (n_nodes: number): Node[] => {
    const root = document.createElement(`section`)
    const nodes: Node[] = [root]
    for (let idx = 1; idx < n_nodes; idx++) {
      const parents = nodes.filter((node) => node instanceof Element)
      const parent = parents[rand_int(parents.length)]
      const kind = rand_int(10)
      const child =
        kind < 6
          ? document.createElement(tags[rand_int(tags.length)])
          : kind < 9
            ? document.createTextNode(`text ${idx}`)
            : document.createComment(`comment ${idx}`)
      const siblings = parent.childNodes
      const next_sibling = siblings[rand_int(siblings.length + 1)]
      if (next_sibling) next_sibling.before(child)
      else parent.append(child)
      nodes.push(child)
    }
    return nodes
  }
  const attached = make_tree(120)
  document.body.append(attached[0])
  const host = attached.find((node) => node instanceof HTMLDivElement)
  if (!host) throw new Error(`random tree has no div to host a shadow root`)
  const shadow_root = host.attachShadow({ mode: `open` })
  const shadow_child = document.createElement(`span`)
  shadow_root.append(shadow_child, document.createTextNode(`shadow text`))
  const detached = make_tree(30) // never appended: a separate tree
  const nodes = [document, document.body, ...attached, shadow_root, shadow_child]
  nodes.push(...detached)

  let n_following = 0
  for (const node_a of nodes) {
    for (const node_b of nodes) {
      const [tree_a, tree_b] = [node_a, node_b] as unknown as HappyDomNode[]
      const expected = walk_is_following(tree_a, tree_b)
      if (NodeUtility.isFollowing(tree_a, tree_b) !== expected) {
        throw new Error(
          `isFollowing mismatch for ${node_a.nodeName} vs ${node_b.nodeName}`,
        )
      }
      n_following += Number(expected)
    }
  }
  // both outcomes are exercised (roughly half the same-tree pairs follow each other)
  expect(n_following).toBeGreaterThan(nodes.length)
  expect(n_following).toBeLessThan(nodes.length ** 2 / 2)
})
