## Basic ring and connectivity information.
##
## Ring membership is derived from bridge detection: a bond lies in a ring
## iff it is not a bridge. Full SSSR perception is not implemented yet.

import std/[algorithm, sequtils]
import ./mol

func connectedComponents*(m: Mol; bondMask: seq[bool] = @[]): seq[int] =
  ## Component id for every atom (ids numbered 0.. in order of first atom).
  ## With a non-empty `bondMask`, only bonds marked `true` connect atoms
  ## (e.g. `ringBonds(m)` gives ring systems).
  result = newSeqWith(m.numAtoms, -1)
  var nextId = 0
  var stack: seq[int]
  for start in 0..<m.numAtoms:
    if result[start] >= 0: continue
    result[start] = nextId
    stack.add start
    while stack.len > 0:
      let a = stack.pop
      for (nb, bi) in m.adj[a]:
        if result[nb] < 0 and (bondMask.len == 0 or bondMask[bi]):
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

func ringBondIndices*(m: Mol; ring: seq[int]): seq[int] =
  ## Bond indices of a ring given as atoms in ring order (as from `sssr`).
  for k in 0..<ring.len:
    result.add m.bondBetween(ring[k], ring[(k + 1) mod ring.len])

# ---------------------------------------------------------------------------
# Smallest set of smallest rings (SSSR)

type RingGraph = object
  ## The ring-bond subgraph with local atom indices.
  atoms: seq[int]                 ## local -> global atom index
  adj: seq[seq[(int, int)]]       ## local neighbour, global bond index
  bondCol: seq[int]               ## global bond -> column in cycle bit vectors (-1 if not a ring bond)
  nCols: int

func buildRingGraph(m: Mol): RingGraph =
  let inRing = ringBonds(m)
  var local = newSeqWith(m.numAtoms, -1)
  result.bondCol = newSeqWith(m.numBonds, -1)
  for bi, bd in m.bonds:
    if not inRing[bi]: continue
    result.bondCol[bi] = result.nCols
    inc result.nCols
    for a in [bd.a, bd.b]:
      if local[a] < 0:
        local[a] = result.atoms.len
        result.atoms.add a
        result.adj.add newSeq[(int, int)]()
    result.adj[local[bd.a]].add (local[bd.b], bi)
    result.adj[local[bd.b]].add (local[bd.a], bi)

func bfsParents(g: RingGraph; src, skipBond: int): seq[(int, int)] =
  ## BFS tree from `src`: (parent, bond) per local atom; (-1, -1) for the
  ## root, (-2, -2) for unreached atoms.
  result = newSeqWith(g.atoms.len, (-2, -2))
  result[src] = (-1, -1)
  var queue = @[src]
  var head = 0
  while head < queue.len:
    let a = queue[head]
    inc head
    for (nb, bi) in g.adj[a]:
      if bi == skipBond or result[nb][0] != -2: continue
      result[nb] = (a, bi)
      queue.add nb

func pathToRoot(parents: seq[(int, int)]; a: int): (seq[int], seq[int]) =
  ## (atoms, bonds) on the tree path from `a` up to the root (atoms include both ends).
  var atoms = @[a]
  var bonds: seq[int]
  var cur = a
  while parents[cur][0] >= 0:
    bonds.add parents[cur][1]
    cur = parents[cur][0]
    atoms.add cur
  (atoms, bonds)

func orderRingAtoms(m: Mol; bonds: seq[int]): seq[int] =
  ## Atom indices of a simple cycle, given its bonds, in walking order.
  var remaining = bonds
  let first = m.bonds[remaining.pop]
  result = @[first.a, first.b]
  while remaining.len > 0:
    let last = result[^1]
    var found = -1
    for k, bi in remaining:
      let bd = m.bonds[bi]
      if bd.a == last or bd.b == last:
        found = k
        break
    let bd = m.bonds[remaining[found]]
    remaining.delete found
    let next = bd.otherAtom(last)
    if remaining.len > 0: result.add next   # last bond closes back to the start

proc sssr*(m: Mol): seq[seq[int]] =
  ## Smallest set of smallest rings: `numRings(m)` independent rings of
  ## minimal total size, each as atom indices in ring order. For symmetric
  ## systems the choice between equal-sized rings is not unique.
  let need = numRings(m)
  if need == 0: return
  let g = buildRingGraph(m)
  var candidates: seq[seq[int]]   # each a list of global bond indices
  # Candidates 1: the shortest cycle through each ring bond.
  for li in 0..<g.atoms.len:
    for (lj, bi) in g.adj[li]:
      if lj < li: continue
      let parents = bfsParents(g, li, bi)
      if parents[lj][0] == -2: continue
      var cyc = pathToRoot(parents, lj)[1]
      cyc.add bi
      candidates.add cyc
  let words = (g.nCols + 63) div 64

  proc select(cands: seq[seq[int]]): seq[seq[int]] =
    ## Greedy GF(2) elimination: shortest linearly independent cycles.
    var sorted = cands
    sorted.sort(proc (x, y: seq[int]): int = cmp(x.len, y.len))
    var basis: seq[(int, seq[uint64])]   # (pivot column, reduced vector)
    for c in sorted:
      var v = newSeq[uint64](words)
      for bi in c:
        let col = g.bondCol[bi]
        v[col div 64] = v[col div 64] xor (1'u64 shl (col mod 64))
      for (piv, bv) in basis:
        if (v[piv div 64] shr (piv mod 64) and 1) == 1:
          for w in 0..<words: v[w] = v[w] xor bv[w]
      var piv = -1
      for w in 0..<words:
        if v[w] != 0:
          for bit in 0..63:
            if (v[w] shr bit and 1) == 1:
              piv = w * 64 + bit
              break
          break
      if piv < 0: continue
      basis.add (piv, v)
      result.add c
      if result.len == need: break

  var chosen = select(candidates)
  if chosen.len < need:
    # Candidates 2 (Horton): P(x,u) + (u,v) + P(x,v) for every ring atom x
    # and ring bond (u,v) whose two tree paths only meet at x.
    for x in 0..<g.atoms.len:
      let parents = bfsParents(g, x, -1)
      for u in 0..<g.atoms.len:
        if parents[u][0] == -2: continue   # other ring system
        for (v, bi) in g.adj[u]:
          if v < u or parents[u][1] == bi or parents[v][1] == bi: continue
          let (atomsU, bondsU) = pathToRoot(parents, u)
          let (atomsV, bondsV) = pathToRoot(parents, v)
          var disjoint = true
          for a in atomsU:
            if a != x and a in atomsV:
              disjoint = false
              break
          if disjoint:
            candidates.add bondsU & bondsV & @[bi]
    chosen = select(candidates)
  for c in chosen:
    result.add orderRingAtoms(m, c)
