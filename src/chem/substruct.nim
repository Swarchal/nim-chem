## Substructure search: matching a `Query` (from SMARTS or a molecule)
## against a molecule by backtracking subgraph isomorphism.
##
## Matches are non-induced (extra bonds between matched atoms are allowed)
## and injective. A match is a sequence giving, for each query atom, the
## index of the molecule atom it maps to. The target should be sanitized
## (aromaticity perceived), as `parseSmiles` does by default.

import std/[algorithm, sets, sequtils]
import ./[mol, valence, rings, aromaticity, query, smarts]
import ./smiles/chirality

type
  Matcher = object
    q: ptr Query
    order: seq[int]      ## query atoms in search order
    parent: seq[int]     ## per query atom: an earlier neighbour, or -1
    map: seq[int]        ## query atom -> target atom (-1 unmapped)
    used: seq[bool]      ## target atom already mapped
    useChirality: bool

  Target = object
    ## The molecule being searched (borrowed for the duration of one
    ## search) and the per-atom properties the query needs, computed once.
    m: ptr Mol
    hCount: seq[int]        ## hydrogens, including H atoms (`H`)
    hCountNoAtoms: seq[int] ## hydrogens not stored as atoms (`h`, `X`)
    valence, ringCount, minRing, ringBondCount: seq[int]
    ringBond: seq[bool]
    recursive: seq[seq[int8]]   ## per recursive query, per atom: -1 unknown, 0/1
    recMatchers: seq[Matcher]   ## per recursive query (empty `order` until used)
    useChirality: bool          ## also applies inside recursive queries

proc initTarget(m: Mol; needs: set[QueryNeeds]; numRecursive: int): Target =
  let n = m.numAtoms
  result.m = addr m
  if qnHydrogens in needs:
    result.hCount = newSeq[int](n)
    result.hCountNoAtoms = newSeq[int](n)
    for i in 0..<n:
      result.hCount[i] = totalHCount(m, i, includeNeighbors = true)
      result.hCountNoAtoms[i] = totalHCount(m, i)
  if needs * {qnRingBonds, qnRingSets} != {}:
    result.ringBond = ringBonds(m)
    result.ringBondCount = newSeq[int](n)
    for bi, bd in m.bonds:
      if result.ringBond[bi]:
        inc result.ringBondCount[bd.a]
        inc result.ringBondCount[bd.b]
  if qnRingSets in needs:
    result.ringCount = newSeq[int](n)
    result.minRing = newSeq[int](n)
    for ring in symmetrizedSssr(m):
      for a in ring:
        inc result.ringCount[a]
        if result.minRing[a] == 0 or ring.len < result.minRing[a]:
          result.minRing[a] = ring.len
  if qnValence in needs:
    # aromatic bonds have no integer order; count them in a Kekulé form
    var k = m
    try: kekulize(k)
    except ChemError: k = m
    result.valence = newSeq[int](n)
    for i in 0..<n: result.valence[i] = totalValence(k, i)
  result.recursive = newSeqWith(numRecursive, newSeqWith(n, -1'i8))
  result.recMatchers = newSeq[Matcher](numRecursive)

proc matchesAt(t: var Target; p: AtomPrim; atom: int): bool

proc atomPrimMatches(t: var Target; p: AtomPrim; i: int): bool =
  let a = t.m.atoms[i]
  case p.kind
  of apAny: true
  of apAtomicNum: a.element == p.value
  of apAromatic: a.aromatic
  of apAliphatic: not a.aromatic
  of apIsotope: a.isotope == p.value
  of apCharge: a.charge == p.value
  of apTotalH: t.hCount[i] == p.value
  of apImplicitH: t.hCountNoAtoms[i] == p.value
  of apDegree: t.m.adj[i].len == p.value
  of apConnectivity: t.m.adj[i].len + t.hCountNoAtoms[i] == p.value
  of apValence: t.valence[i] == p.value
  of apRingCount: t.ringCount[i] == p.value
  of apRingSize: t.minRing[i] == p.value
  of apRingBonds: t.ringBondCount[i] == p.value
  of apInRing: t.ringBondCount[i] > 0
  of apRecursive:
    if t.recursive[p.subId][i] < 0:
      t.recursive[p.subId][i] = int8(matchesAt(t, p, i))
    t.recursive[p.subId][i] == 1

proc atomMatches(t: var Target; e: Expr[AtomPrim]; i: int): bool =
  case e.op
  of eoPrim: atomPrimMatches(t, e.prim, i)
  of eoNot: not atomMatches(t, e.arg, i)
  of eoAnd: e.args.allIt(atomMatches(t, it, i))
  of eoOr: e.args.anyIt(atomMatches(t, it, i))

func bondMatches(t: Target; e: Expr[BondPrim]; bi: int): bool =
  case e.op
  of eoPrim:
    let order = t.m.bonds[bi].order
    case e.prim
    of bpDefault: order in {boSingle, boAromatic}
    of bpSingle: order == boSingle
    of bpDouble: order == boDouble
    of bpTriple: order == boTriple
    of bpQuadruple: order == boQuadruple
    of bpAromatic: order == boAromatic
    of bpAny: true
    of bpRing: t.ringBond[bi]
  of eoNot: not bondMatches(t, e.arg, bi)
  of eoAnd: e.args.allIt(bondMatches(t, it, bi))
  of eoOr: e.args.anyIt(bondMatches(t, it, bi))

proc initMatcher(q: ptr Query; nTarget: int; useChirality: bool): Matcher =
  ## Search order: breadth-first from atom 0 through each component, so
  ## every atom but a component's first has an already-placed neighbour.
  ## `q` must outlive the matcher.
  result.q = q
  result.parent = newSeqWith(q[].numAtoms, -1)
  result.map = newSeqWith(q[].numAtoms, -1)
  result.used = newSeq[bool](nTarget)
  result.useChirality = useChirality
  var seen = newSeq[bool](q[].numAtoms)
  for root in 0..<q[].numAtoms:
    if seen[root]: continue
    seen[root] = true
    var head = result.order.len
    result.order.add root
    while head < result.order.len:
      let a = result.order[head]
      inc head
      for (nb, _) in q.adj[a]:
        if not seen[nb]:
          seen[nb] = true
          result.parent[nb] = a
          result.order.add nb

proc stereoMatches(mt: Matcher; t: Target): bool =
  ## Every query atom with `@`/`@@` maps to a target atom with the same
  ## spatial arrangement of the matched neighbours, and every query double
  ## bond with cis/trans to a target bond with the same arrangement.
  for qa in 0..<mt.q[].numAtoms:
    let qatom = mt.q.atoms[qa]
    if qatom.chiral == ctNone: continue
    let ta = mt.map[qa]
    if t.m.atoms[ta].chiral == ctNone: return false
    var written = qatom.nbrOrder.mapIt(if it < 0: -1 else: mt.map[it])
    let reference = chiralReference(t.m[], ta)
    var extra: seq[int]   # target positions the query doesn't mention
    for x in reference:
      if x notin written: extra.add x
    if extra.len == 1:
      let k = written.find(-1)
      if k >= 0: written[k] = extra[0]   # query H matched by an H atom
      else: written.add extra[0]         # unmentioned neighbour goes last (as RDKit)
    elif extra.len > 1:
      continue   # too few query neighbours to define a configuration
    if written.len != reference.len or written.anyIt(it notin reference):
      continue
    let tag = reorderChirality(t.m[], ta, qatom.chiral, written, reference,
                               qatom.ringDigits)
    if tag != t.m.atoms[ta].chiral: return false
  for qbond in mt.q.bonds:
    if qbond.stereo == bsNone: continue
    let (ta, tb) = (mt.map[qbond.a], mt.map[qbond.b])
    let tbond = t.m.bonds[t.m[].bondBetween(ta, tb)]
    if tbond.stereo == bsNone: return false
    var refs = (mt.map[qbond.stereoAtoms[0]], mt.map[qbond.stereoAtoms[1]])
    if tbond.a != ta: refs = (refs[1], refs[0])
    # the other substituent on an end is on the opposite side
    var cis = tbond.stereo == bsCis
    if refs[0] != tbond.stereoAtoms[0]: cis = not cis
    if refs[1] != tbond.stereoAtoms[1]: cis = not cis
    if cis != (qbond.stereo == bsCis): return false
  true

proc search(mt: var Matcher; t: var Target; depth: int;
            onMatch: proc (map: seq[int]): bool): bool =
  ## Extends the partial match from `order[depth]`. `onMatch` is called
  ## for each complete match and returns false to stop the search; the
  ## result is false once stopped.
  if depth == mt.order.len:
    if mt.useChirality and not stereoMatches(mt, t): return true
    return onMatch(mt.map)
  let qa = mt.order[depth]
  let par = mt.parent[qa]

  template tryAtom(ta, parentBond: int) =
    ## `parentBond`: the bond from the parent's image to `ta`, if known
    if not mt.used[ta] and atomMatches(t, mt.q.atoms[qa].expr, ta):
      var ok = true
      for (qn, qb) in mt.q.adj[qa]:
        let tn = mt.map[qn]
        if tn < 0: continue
        let tb = if qn == par: parentBond else: t.m[].bondBetween(ta, tn)
        if tb < 0 or not bondMatches(t, mt.q.bonds[qb].expr, tb):
          ok = false
          break
      if ok:
        mt.map[qa] = ta
        mt.used[ta] = true
        let more = search(mt, t, depth + 1, onMatch)
        mt.map[qa] = -1
        mt.used[ta] = false
        if not more: return false

  if par < 0:
    for ta in 0..<t.m[].numAtoms: tryAtom(ta, -1)
  else:
    for (ta, tb) in t.m.adj[mt.map[par]]: tryAtom(ta, tb)
  true

proc matchesAt(t: var Target; p: AtomPrim; atom: int): bool =
  ## Whether recursive query `p.sub` matches with its first atom on `atom`.
  if t.recMatchers[p.subId].order.len == 0:
    t.recMatchers[p.subId] = initMatcher(addr p.sub[], t.m[].numAtoms, t.useChirality)
  if not atomMatches(t, p.sub.atoms[0].expr, atom): return false
  # the matcher is moved out while searching: the search may evaluate
  # other recursive queries, which needs `t` (it can't reach this one)
  var mt = move t.recMatchers[p.subId]
  mt.map[0] = atom
  mt.used[atom] = true
  var found = false
  discard search(mt, t, 1, proc (map: seq[int]): bool =
    found = true
    false)
  mt.map[0] = -1
  mt.used[atom] = false
  t.recMatchers[p.subId] = move mt
  found

proc substructMatches*(m: Mol; q: Query; uniquify = true; maxMatches = 1000;
                       useChirality = false): seq[seq[int]] =
  ## All matches of `q` in `m` (up to `maxMatches`), each giving the
  ## molecule atom for every query atom. With `uniquify` (default), matches
  ## covering the same set of atoms are reported once. With `useChirality`,
  ## query `@`/`@@` and `/` `\` cis/trans must agree with the molecule's
  ## stereo, which must be specified there (otherwise stereo is ignored).
  if q.numAtoms == 0 or maxMatches <= 0: return
  var t = initTarget(m, q.needs, q.numRecursive)
  t.useChirality = useChirality
  var mt = initMatcher(addr q, m.numAtoms, useChirality)
  var seen: HashSet[seq[int]]
  var found: seq[seq[int]]
  var key = newSeq[int](q.numAtoms)   # reused: copied only when new
  discard search(mt, t, 0, proc (map: seq[int]): bool =
    if uniquify:
      key[0..^1] = map
      key.sort()
      if key in seen: return true
      seen.incl key
    found.add map
    found.len < maxMatches)
  found

proc substructMatch*(m: Mol; q: Query; useChirality = false): seq[int] =
  ## The first match of `q` in `m`, or an empty sequence.
  let ms = substructMatches(m, q, uniquify = false, maxMatches = 1,
                            useChirality = useChirality)
  if ms.len > 0: ms[0] else: @[]

proc hasSubstructMatch*(m: Mol; q: Query; useChirality = false): bool =
  ## Whether `q` occurs in `m`.
  substructMatch(m, q, useChirality).len > 0

func asQuery(q: Mol): Query = toQuery(q)
proc asQuery(q: string): Query = parseSmarts(q)

proc substructMatches*(m: Mol; q: Mol | string; uniquify = true;
                       maxMatches = 1000; useChirality = false): seq[seq[int]] =
  ## Convenience: a molecule query (see `query.toQuery` for how atoms
  ## compare) or a SMARTS string (see `smarts.parseSmarts`).
  substructMatches(m, asQuery(q), uniquify, maxMatches, useChirality)

proc substructMatch*(m: Mol; q: Mol | string; useChirality = false): seq[int] =
  substructMatch(m, asQuery(q), useChirality)

proc hasSubstructMatch*(m: Mol; q: Mol | string; useChirality = false): bool =
  hasSubstructMatch(m, asQuery(q), useChirality)
