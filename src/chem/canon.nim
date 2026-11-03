## Canonical atom ranking by iterative partition refinement
## (Weininger CANON / Morgan style), with tie-breaking for symmetric atoms.
##
## Stereochemistry is not part of the invariants (v1): ranks depend only on
## the constitution, so stereoisomers share an atom ordering and differ only
## in their written `@`/`@@` and `/`/`\` markers.

import std/[algorithm, sequtils]
import ./[mol, valence, rings]

func cmpSeq(a, b: seq[int]): int =
  for i in 0..<min(a.len, b.len):
    if a[i] != b[i]: return cmp(a[i], b[i])
  cmp(a.len, b.len)

func denseRanks(keys: seq[seq[int]]): (seq[int], int) =
  ## Dense ranks (0..k-1) of keys, equal keys sharing a rank, and k.
  var idx = toSeq(0..<keys.len)
  idx.sort(proc (x, y: int): int = cmpSeq(keys[x], keys[y]))
  var ranks = newSeq[int](keys.len)
  var r = 0
  for k, i in idx:
    if k > 0 and cmpSeq(keys[idx[k-1]], keys[i]) != 0: inc r
    ranks[i] = r
  (ranks, if keys.len == 0: 0 else: r + 1)

proc refine(m: Mol; ranks: var seq[int]; nClasses: var int) =
  ## Splits rank classes by sorted neighbour (rank, bond order) codes until
  ## the partition is stable.
  var keys = newSeq[seq[int]](m.numAtoms)
  while true:
    for i in 0..<m.numAtoms:
      var key = @[ranks[i]]
      for (n, bi) in m.adj[i]:
        key.add ranks[n] * 8 + ord(m.bonds[bi].order)
      sort(key.toOpenArray(1, key.high))
      keys[i] = key
    let (newRanks, nd) = denseRanks(keys)
    ranks = newRanks
    if nd == nClasses: break
    nClasses = nd

proc canonicalRanks*(m: Mol): seq[int] =
  ## A canonical, unique rank (0..n-1) for every atom. Lower ranks are
  ## written first by the canonical SMILES writer.
  let n = m.numAtoms
  if n == 0: return @[]
  let inRing = ringBonds(m)
  var keys = newSeq[seq[int]](n)
  for i, a in m.atoms:
    var ringDeg = 0
    for (_, bi) in m.adj[i]:
      if inRing[bi]: inc ringDeg
    keys[i] = @[m.degree(i), a.element, a.isotope, a.charge,
                totalHCount(m, i), ord(a.aromatic), ringDeg, a.atomMap]
  var nClasses: int
  (result, nClasses) = denseRanks(keys)
  refine(m, result, nClasses)
  # Break ties: split the lowest tied class by keeping its first atom and
  # moving the rest (and every higher class) up one rank.
  while nClasses < n:
    var seen = newSeq[bool](n)
    var tied = -1
    for r in result:
      if seen[r]:
        if tied < 0 or r < tied: tied = r
      seen[r] = true
    let chosen = result.find(tied)
    for i in 0..<n:
      if result[i] > tied or (result[i] == tied and i != chosen): inc result[i]
    inc nClasses
    refine(m, result, nClasses)
