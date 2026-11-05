## Canonical atom ranking by iterative partition refinement
## (Weininger CANON / Morgan style), with tie-breaking for symmetric atoms.
##
## Stereochemistry is folded in as invariants (after Schneider, Sayle and
## Landrum, J. Chem. Inf. Model. 2015): once the neighbours of a stereocentre
## (or the substituents of a stereo double bond) have distinct ranks, its
## configuration relative to those ranks becomes a label that can split
## otherwise-equivalent atoms. This runs after every refinement and every
## tie-break, so e.g. the two centres of a meso compound are ordered by
## their (opposite) labels instead of by input order.

import std/[algorithm, sequtils]
import ./[mol, valence, rings]

const maxSearchRefinements = 5000
  ## cap on tie-break candidates tried while searching stereo tie-breaks

func cmpSeq(a, b: seq[int]): int =
  for i in 0 ..< min(a.len, b.len):
    if a[i] != b[i]:
      return cmp(a[i], b[i])
  cmp(a.len, b.len)

func denseRanks(keys: seq[seq[int]]): (seq[int], int) =
  ## Dense ranks (0..k-1) of keys, equal keys sharing a rank, and k.
  var idx = toSeq(0 ..< keys.len)
  idx.sort(
    proc(x, y: int): int =
      cmpSeq(keys[x], keys[y])
  )
  var ranks = newSeq[int](keys.len)
  var r = 0
  for k, i in idx:
    if k > 0 and cmpSeq(keys[idx[k - 1]], keys[i]) != 0:
      inc r
    ranks[i] = r
  (
    ranks,
    if keys.len == 0:
      0
    else:
      r + 1,
  )

proc refine(m: Mol, ranks: var seq[int], nClasses: var int) =
  ## Splits rank classes by sorted neighbour (rank, bond order) codes until
  ## the partition is stable.
  var keys = newSeq[seq[int]](m.numAtoms)
  while true:
    for i in 0 ..< m.numAtoms:
      var key = @[ranks[i]]
      for (n, bi) in m.adj[i]:
        key.add ranks[n] * 8 + ord(m.bonds[bi].order)
      sort(key.toOpenArray(1, key.high))
      keys[i] = key
    let (newRanks, nd) = denseRanks(keys)
    ranks = newRanks
    if nd == nClasses:
      break
    nClasses = nd

func tetrahedralLabel(m: Mol, i: int, ranks: seq[int]): int =
  ## 1 or 2 for a chiral atom whose stereo positions all have distinct ranks
  ## (the implicit H / lone pair ranks lowest), else 0.
  let tag = m.atoms[i].chiral
  if tag == ctNone:
    return 0
  let reference = chiralReference(m, i)
  if reference.len < 3:
    return 0
  var byRank = newSeq[(int, int)](reference.len) # (rank, atom)
  for k, x in reference:
    byRank[k] = ((if x < 0: -1 else: ranks[x]), x)
  byRank.sort()
  for k in 1 ..< byRank.len:
    if byRank[k][0] == byRank[k - 1][0]:
      return 0
  (ord(tag == ctCW) xor permutationParity(reference, byRank.mapIt(it[1]))) + 1

func substituents*(m: Mol, endAtom, other: int): seq[int] =
  ## Neighbours of `endAtom` other than `other` (its double-bond partner).
  for nb in m.neighbors(endAtom):
    if nb != other:
      result.add nb

func doubleBondLabel(m: Mol, bi: int, ranks: seq[int]): int =
  ## 1 (cis) or 2 (trans) relative to the lowest-ranked substituent at each
  ## end, once each end's substituents have distinct ranks, else 0.
  let bd = m.bonds[bi]
  if bd.stereo == bsNone:
    return 0
  var cis = bd.stereo == bsCis
  for (endAtom, other, refAtom) in [
    (bd.a, bd.b, bd.stereoAtoms[0]), (bd.b, bd.a, bd.stereoAtoms[1])
  ]:
    let subs = substituents(m, endAtom, other)
    if subs.len == 2:
      if ranks[subs[0]] == ranks[subs[1]]:
        return 0
      let lowest =
        if ranks[subs[0]] < ranks[subs[1]]:
          subs[0]
        else:
          subs[1]
      if lowest != refAtom:
        cis = not cis
    elif subs.len != 1:
      return 0
  if cis: 1 else: 2

func stereoKeys(m: Mol, ranks: seq[int], stereoBonds: seq[int]): seq[seq[int]] =
  ## (rank, tetrahedral label, double-bond label) per atom.
  result = newSeq[seq[int]](m.numAtoms)
  for i in 0 ..< m.numAtoms:
    result[i] = @[ranks[i], tetrahedralLabel(m, i, ranks), 0]
  for bi in stereoBonds:
    let label = doubleBondLabel(m, bi, ranks)
    result[m.bonds[bi].a][2] = label
    result[m.bonds[bi].b][2] = label

func stereoBondList(m: Mol): seq[int] =
  for bi, bd in m.bonds:
    if bd.stereo != bsNone:
      result.add bi

func hasStereo*(m: Mol): bool =
  ## Whether any atom has a chirality tag or any bond cis/trans stereo.
  m.atoms.anyIt(it.chiral != ctNone) or m.bonds.anyIt(it.stereo != bsNone)

func allLabelled(m: Mol, ranks: seq[int], stereoBonds: seq[int]): bool =
  ## Every stereo element has a label under `ranks`.
  for i in 0 ..< m.numAtoms:
    if m.atoms[i].chiral != ctNone and tetrahedralLabel(m, i, ranks) == 0:
      return false
  stereoBonds.allIt(doubleBondLabel(m, it, ranks) != 0)

proc refineWithStereo(
    m: Mol,
    ranks: var seq[int],
    nClasses: var int,
    useStereo: bool,
    stereoBonds: seq[int],
) =
  ## `refine`, then (with `useStereo`) repeatedly split classes by stereo
  ## labels that the current ranks make computable.
  refine(m, ranks, nClasses)
  if not useStereo:
    return
  while nClasses < m.numAtoms:
    let (newRanks, nd) = denseRanks(stereoKeys(m, ranks, stereoBonds))
    if nd == nClasses:
      break
    ranks = newRanks
    nClasses = nd
    refine(m, ranks, nClasses)

proc initialRanks(m: Mol): (seq[int], int) =
  let n = m.numAtoms
  let inRing = ringBonds(m)
  var keys = newSeq[seq[int]](n)
  for i, a in m.atoms:
    var ringDeg = 0
    for (_, bi) in m.adj[i]:
      if inRing[bi]:
        inc ringDeg
    keys[i] = @[
      m.degree(i),
      a.element,
      a.isotope,
      a.charge,
      totalHCount(m, i),
      ord(a.aromatic),
      ringDeg,
      a.atomMap,
    ]
  denseRanks(keys)

proc symmetryClasses*(m: Mol, useStereo = true): seq[int] =
  ## Symmetry classes: atoms with equal values are indistinguishable by
  ## constitution (and, with `useStereo`, by computable stereo labels).
  ## No tie-breaking, so values are not unique.
  if m.numAtoms == 0:
    return @[]
  var nClasses: int
  (result, nClasses) = initialRanks(m)
  let stereo = useStereo and hasStereo(m)
  refineWithStereo(
    m,
    result,
    nClasses,
    stereo,
    if stereo:
      stereoBondList(m)
    else:
      @[],
  )

proc canonicalRanks*(m: Mol, useStereo = true): seq[int] =
  ## A canonical, unique rank (0..n-1) for every atom. Lower ranks are
  ## written first by the canonical SMILES writer.
  let n = m.numAtoms
  if n == 0:
    return @[]
  var nClasses: int
  (result, nClasses) = initialRanks(m)
  let stereo = useStereo and hasStereo(m)
  let stereoBonds =
    if stereo:
      stereoBondList(m)
    else:
      @[]
  refineWithStereo(m, result, nClasses, stereo, stereoBonds)
  if nClasses == n:
    return

  proc breakTie(ranks: seq[int], nClasses, tied, chosen: int): (seq[int], int) =
    ## Splits class `tied` by promoting atom `chosen` (the rest, and every
    ## higher class, move up one rank), then refines.
    var r = ranks
    for i in 0 ..< n:
      if r[i] > tied or (r[i] == tied and i != chosen):
        inc r[i]
    var nc = nClasses + 1
    refineWithStereo(m, r, nc, stereo, stereoBonds)
    (r, nc)

  proc lowestTied(ranks: seq[int]): int =
    var seen = newSeq[bool](n)
    result = -1
    for r in ranks:
      if seen[r] and (result < 0 or r < result):
        result = r
      seen[r] = true

  proc signature(ranks: seq[int]): seq[int] =
    ## stereo labels in rank order: compares the outcomes of tie-breaks
    var keys = stereoKeys(m, ranks, stereoBonds)
    keys.sort(cmpSeq)
    for k in keys:
      result.add k[1]
      result.add k[2]

  proc encoding(ranks: seq[int]): seq[int] =
    ## The molecule written in terms of complete ranks (bonds, then stereo
    ## labels): equal encodings mean the same labelled structure.
    var bonds = m.bonds.mapIt(
      (min(ranks[it.a], ranks[it.b]), max(ranks[it.a], ranks[it.b]), ord(it.order))
    )
    bonds.sort()
    for (x, y, o) in bonds:
      result.add [x, y, o]
    result.add signature(ranks)

  # Without stereo, or once every stereo element has a label, the choice of
  # tied atom can't change any label (tie-breaks keep the order of distinct
  # ranks), so the first is taken. Before that it can (e.g. the two ring
  # directions of cis-1,4-dimethylcyclohexane flip both centres' labels):
  # candidates are compared by their stereo labels, and if several tie,
  # each is completed and the smallest final encoding wins (a search over
  # an atom-order-independent set of choices, so the result is canonical).
  # `budget` caps the number of refinements spent searching; beyond it
  # the first candidate is taken.
  var budget = maxSearchRefinements

  proc complete(ranks0: seq[int], nc0: int): seq[int] =
    var (ranks, nc) = (ranks0, nc0)
    while nc < n:
      let tied = lowestTied(ranks)
      if not stereo or budget <= 0 or allLabelled(m, ranks, stereoBonds):
        (ranks, nc) = breakTie(ranks, nc, tied, ranks.find(tied))
        continue
      var best: seq[(seq[int], int)]
      var bestSig: seq[int]
      for i in 0 ..< n:
        if ranks[i] != tied:
          continue
        dec budget
        let candidate = breakTie(ranks, nc, tied, i)
        let sig = signature(candidate[0])
        let c =
          if best.len == 0:
            -1
          else:
            cmpSeq(sig, bestSig)
        if c < 0:
          best = @[candidate]
          bestSig = sig
        elif c == 0:
          best.add candidate
      if best.len == 1:
        (ranks, nc) = best[0]
        continue
      var bestRanks: seq[int]
      var bestCode: seq[int]
      for (r, k) in best:
        let full = complete(r, k)
        let code = encoding(full)
        if bestRanks.len == 0 or cmpSeq(code, bestCode) < 0:
          bestRanks = full
          bestCode = code
      return bestRanks
    ranks

  result = complete(result, nClasses)
