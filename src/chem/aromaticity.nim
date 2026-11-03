## Kekulization and aromaticity perception.
##
## `kekulize` turns aromatic bonds into an explicit single/double pattern
## (a perfect matching found with Edmonds' blossom algorithm), and
## `setAromaticity` re-derives aromatic flags from a Kekulé structure using a
## Hückel 4n+2 rule over SSSR rings and fused ring combinations, with
## RDKit-style electron counting. `sanitize` runs both plus valence checks,
## so equivalent aromatic and Kekulé inputs end up identical.
##
## Hydrogen counts of atoms that are (or become) aromatic are frozen into
## `explicitH` (with `noImplicit = true`), so they no longer depend on the
## bond pattern; this is what keeps e.g. pyrrole's N-H when its bonds change.

import std/[sequtils, sets, algorithm, tables]
import ./[mol, valence, rings]

type KekulizeError* = object of ChemError

const maxFusedRings = 6
  ## largest number of fused rings considered together as one pi system

func electronegativity(z: int): float =
  ## Pauling electronegativity for elements that occur in aromatic rings
  ## or their exocyclic double bonds (0 for others).
  case z
  of 1: 2.20
  of 5: 2.04
  of 6: 2.55
  of 7: 3.04
  of 8: 3.44
  of 9: 3.98
  of 14: 1.90
  of 15: 2.19
  of 16: 2.58
  of 17: 3.16
  of 33: 2.18
  of 34: 2.55
  of 35: 2.96
  of 52: 2.10
  of 53: 2.66
  else: 0.0

proc freezeHydrogens(m: var Mol; i: int) =
  if not m.atoms[i].noImplicit:
    m.atoms[i].explicitH = totalHCount(m, i)
    m.atoms[i].noImplicit = true

proc maxMatching(g: seq[seq[int]]): seq[int] =
  ## Maximum matching in a general graph (Edmonds' blossom algorithm).
  ## Returns the partner of each vertex, or -1.
  let n = g.len
  var match = newSeqWith(n, -1)
  for v in 0..<n:   # greedy start
    if match[v] == -1:
      for u in g[v]:
        if match[u] == -1:
          match[u] = v
          match[v] = u
          break
  var parent = newSeq[int](n)
  var base = newSeq[int](n)
  var used = newSeq[bool](n)
  var blossom = newSeq[bool](n)

  proc lca(a0, b0: int): int =
    var seen = newSeq[bool](n)
    var a = a0
    while true:
      a = base[a]
      seen[a] = true
      if match[a] == -1: break
      a = parent[match[a]]
    var b = b0
    while true:
      b = base[b]
      if seen[b]: return b
      b = parent[match[b]]

  proc markPath(v0, b, child0: int) =
    var v = v0
    var child = child0
    while base[v] != b:
      blossom[base[v]] = true
      blossom[base[match[v]]] = true
      parent[v] = child
      child = match[v]
      v = parent[match[v]]

  proc findPath(root: int): int =
    for i in 0..<n:
      used[i] = false
      parent[i] = -1
      base[i] = i
    used[root] = true
    var queue = @[root]
    var head = 0
    while head < queue.len:
      let v = queue[head]
      inc head
      for to in g[v]:
        if base[v] == base[to] or match[v] == to: continue
        if to == root or (match[to] != -1 and parent[match[to]] != -1):
          let cur = lca(v, to)
          for i in 0..<n: blossom[i] = false
          markPath(v, cur, to)
          markPath(to, cur, v)
          for i in 0..<n:
            if blossom[base[i]]:
              base[i] = cur
              if not used[i]:
                used[i] = true
                queue.add i
        elif parent[to] == -1:
          parent[to] = v
          if match[to] == -1: return to
          used[match[to]] = true
          queue.add match[to]
    -1

  for v in 0..<n:
    if match[v] == -1:
      var u = findPath(v)
      while u != -1:
        let pv = parent[u]
        let ppv = match[pv]
        match[u] = pv
        match[pv] = u
        u = ppv
  match

proc kekulize*(m: var Mol) =
  ## Replaces every aromatic bond by a single or double bond so that each
  ## aromatic atom that needs one gets exactly one double bond, then clears
  ## all aromatic flags. Raises `KekulizeError` if that is impossible
  ## (e.g. `c1cccc1`, or pyrrole written without its [nH]).
  for i in 0..<m.numAtoms:
    var touches = m.atoms[i].aromatic
    for (_, bi) in m.adj[i]:
      if m.bonds[bi].order == boAromatic: touches = true
    if touches: freezeHydrogens(m, i)
  # An atom needs a double bond if its valence (aromatic bonds counted as 1)
  # is below its lowest allowed valence.
  var local = newSeqWith(m.numAtoms, -1)
  var candidates: seq[int]
  for i in 0..<m.numAtoms:
    var nArom = 0
    var v = m.atoms[i].explicitH
    for (_, bi) in m.adj[i]:
      if m.bonds[bi].order == boAromatic: inc nArom
      v += bondValence(m.bonds[bi].order)
    if nArom == 0: continue
    let vals = allowedValences(m.atoms[i])
    if vals.len > 0 and v < vals[0]:
      local[i] = candidates.len
      candidates.add i
  var g = newSeq[seq[int]](candidates.len)
  for bd in m.bonds:
    if bd.order == boAromatic and local[bd.a] >= 0 and local[bd.b] >= 0:
      g[local[bd.a]].add local[bd.b]
      g[local[bd.b]].add local[bd.a]
  for nbrs in g.mitems: nbrs.sort()   # result depends only on atom order
  let match = maxMatching(g)
  var unmatched: seq[int]
  for li, partner in match:
    if partner < 0: unmatched.add candidates[li]
  if unmatched.len > 0:
    raise newException(KekulizeError,
      "can't kekulize: no double bond possible for atom(s) " & $unmatched)
  for bd in m.bonds.mitems:
    if bd.order != boAromatic: continue
    let la = local[bd.a]
    bd.order = if la >= 0 and match[la] == local[bd.b]: boDouble else: boSingle
  for a in m.atoms.mitems: a.aromatic = false

func piElectrons(m: Mol; i: int; inRing: seq[bool]): int =
  ## Electrons atom `i` donates to a ring pi system in a Kekulé structure,
  ## or -1 if it cannot be part of an aromatic ring.
  let a = m.atoms[i]
  if a.element == 0: return -1
  var ringDouble, exoDouble = 0
  var exoPartner = -1
  for (nb, bi) in m.adj[i]:
    case m.bonds[bi].order
    of boDouble:
      if inRing[bi]: inc ringDouble
      else:
        inc exoDouble
        exoPartner = nb
    of boSingle: discard
    else: return -1   # triple/quadruple/aromatic/unknown
  let h = totalHCount(m, i)
  if m.degree(i) + h > 3: return -1   # not sp2
  if ringDouble == 1 and exoDouble == 0: return 1
  if ringDouble > 0 or exoDouble > 1: return -1
  if exoDouble == 1:
    # exocyclic C=O, C=N, [N+]=O ...: the double bond's electrons are not
    # shared with the ring (RDKit applies this to C and N, not S or P)
    return (if a.element in [6, 7] and electronegativity(m.atoms[exoPartner].element) >
               electronegativity(a.element): 0 else: -1)
  let total = explicitValence(m, i) + h
  let (z, q) = (a.element, a.charge)
  if (z in [7, 15, 33] and q == 0 and total == 3) or
     (z in [8, 16, 34, 52] and q == 0 and total == 2) or
     (z in [8, 16] and q == 1 and total == 3) or
     (z == 6 and q == -1 and total == 3) or
     (z == 7 and q == -1 and total == 2):
    return 2   # lone pair
  if z == 6 and q == 1 and total == 3:
    return 0   # carbocation: empty p orbital (neutral boron is not
               # treated as aromatic, matching RDKit)
  -1

proc connectedRingSubsets(adj: seq[seq[int]]; members: seq[int];
                          maxSize: int): seq[seq[int]] =
  ## All connected subsets (as sorted ring indices) of size 2..maxSize of a
  ## fused ring system.
  var seen = initHashSet[seq[int]]()
  var frontier: seq[seq[int]]
  for r in members: frontier.add @[r]
  for size in 2..maxSize:
    var next: seq[seq[int]]
    for sub in frontier:
      for r in sub:
        for nb in adj[r]:
          if nb in sub: continue
          var grown = sub & nb
          grown.sort()
          if grown notin seen:
            seen.incl grown
            next.add grown
    result.add next
    frontier = next
    if frontier.len == 0: break

proc setAromaticity*(m: var Mol) =
  ## Marks atoms and bonds of aromatic rings (Hückel 4n+2 over single SSSR
  ## rings and fused combinations of up to `maxFusedRings` rings). The
  ## molecule is kekulized first if it still has aromatic bonds; existing
  ## aromatic flags are replaced.
  kekulize(m)
  let rings = sssr(m)
  if rings.len == 0: return
  let inRing = ringBonds(m)
  var electrons = newSeq[int](m.numAtoms)
  for i in 0..<m.numAtoms: electrons[i] = piElectrons(m, i, inRing)
  # ring bonds, and fusion (shared bond) between rings
  var ringBondList = newSeq[seq[int]](rings.len)
  for r, ring in rings:
    for k in 0..<ring.len:
      ringBondList[r].add m.bondBetween(ring[k], ring[(k + 1) mod ring.len])
  var fused = newSeq[seq[int]](rings.len)
  for r1 in 0..<rings.len:
    for r2 in r1+1..<rings.len:
      if ringBondList[r1].anyIt(it in ringBondList[r2]):
        fused[r1].add r2
        fused[r2].add r1

  var aromaticRing = newSeq[bool](rings.len)
  let bondAtoms = m.bonds.mapIt((it.a, it.b))   # closures can't capture `var m`
  proc check(subset: seq[int]): bool =
    ## Hückel test. For fused combinations only perimeter atoms count:
    ## those on a bond that belongs to just one ring of the combination
    ## (so e.g. pyrene's two central atoms are excluded).
    var atoms = initHashSet[int]()
    for r in subset:
      for a in rings[r]:
        if electrons[a] < 0: return false
    var bondUse = initCountTable[int]()
    for r in subset:
      for bi in ringBondList[r]: bondUse.inc bi
    for bi, uses in bondUse:
      if uses == 1:
        atoms.incl bondAtoms[bi][0]
        atoms.incl bondAtoms[bi][1]
    var total = 0
    for a in atoms: total += electrons[a]
    total mod 4 == 2

  for r in 0..<rings.len:
    if check(@[r]): aromaticRing[r] = true
  # fused systems: connected components of the `fused` graph
  var system = newSeqWith(rings.len, -1)
  for start in 0..<rings.len:
    if system[start] >= 0 or fused[start].len == 0: continue
    var members = @[start]
    system[start] = start
    var k = 0
    while k < members.len:
      for nb in fused[members[k]]:
        if system[nb] < 0:
          system[nb] = start
          members.add nb
      inc k
    # only rings whose atoms can all take part are worth combining
    let usable = members.filterIt(rings[it].allIt(electrons[it] >= 0))
    if usable.len < 2: continue
    var usableAdj = newSeq[seq[int]](rings.len)
    for r in usable:
      usableAdj[r] = fused[r].filterIt(it in usable)
    for subset in connectedRingSubsets(usableAdj, usable,
                                       min(maxFusedRings, usable.len)):
      if subset.allIt(aromaticRing[it]): continue
      if check(subset):
        for r in subset: aromaticRing[r] = true

  var aromAtoms: seq[int]
  for r, ring in rings:
    if aromaticRing[r]: aromAtoms.add ring
  for a in aromAtoms: freezeHydrogens(m, a)   # before bond orders change
  for a in aromAtoms: m.atoms[a].aromatic = true
  for r in 0..<rings.len:
    if aromaticRing[r]:
      for bi in ringBondList[r]: m.bonds[bi].order = boAromatic

proc sanitize*(m: var Mol) =
  ## Kekulizes, checks valences on the Kekulé structure, then perceives
  ## aromaticity. Raises `KekulizeError` or `ValenceError`.
  kekulize(m)
  checkValences(m)
  setAromaticity(m)
