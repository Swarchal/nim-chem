## Basic ring and connectivity information.
##
## Ring membership is derived from bridge detection: a bond lies in a ring
## iff it is not a bridge. Full SSSR perception is not implemented yet.

import std/sequtils
import ./mol

func connectedComponents*(m: Mol): seq[int] =
  ## Component id for every atom (ids numbered 0.. in order of first atom).
  result = newSeqWith(m.numAtoms, -1)
  var nextId = 0
  var stack: seq[int]
  for start in 0..<m.numAtoms:
    if result[start] >= 0: continue
    result[start] = nextId
    stack.add start
    while stack.len > 0:
      let a = stack.pop
      for nb in m.neighbors(a):
        if result[nb] < 0:
          result[nb] = nextId
          stack.add nb
    inc nextId

func componentCount*(comp: openArray[int]): int =
  ## Number of components given the result of `connectedComponents`.
  if comp.len == 0: 0 else: max(comp) + 1

func numComponents*(m: Mol): int = componentCount(connectedComponents(m))

func ringBonds*(m: Mol): seq[bool] =
  ## `true` for every bond that is part of at least one ring.
  ## Iterative Tarjan bridge-finding.
  let n = m.numAtoms
  var disc = newSeqWith(n, -1)
  var low = newSeq[int](n)
  result = newSeqWith(m.numBonds, true)
  var t = 0
  # stack frames: (atom, parent bond, next adjacency index)
  var stack: seq[(int, int, int)]
  for root in 0..<n:
    if disc[root] >= 0: continue
    disc[root] = t; low[root] = t; inc t
    stack.add (root, -1, 0)
    while stack.len > 0:
      let (a, pb, idx) = stack[^1]
      if idx < m.adj[a].len:
        stack[^1][2] = idx + 1
        let (nb, bi) = m.adj[a][idx]
        if bi == pb: continue
        if disc[nb] < 0:
          disc[nb] = t; low[nb] = t; inc t
          stack.add (nb, bi, 0)
        else:
          low[a] = min(low[a], disc[nb])
      else:
        discard stack.pop
        if stack.len > 0:
          let p = stack[^1][0]
          low[p] = min(low[p], low[a])
          if low[a] > disc[p]:
            result[pb] = false   # bridge

func ringAtoms*(m: Mol): seq[bool] =
  ## `true` for atoms that belong to at least one ring.
  result = newSeq[bool](m.numAtoms)
  let rb = ringBonds(m)
  for i, bd in m.bonds:
    if rb[i]:
      result[bd.a] = true
      result[bd.b] = true

func numRings*(m: Mol): int =
  ## Number of rings in the smallest set of smallest rings, computed as the
  ## cyclomatic number E - V + C.
  m.numBonds - m.numAtoms + m.numComponents
