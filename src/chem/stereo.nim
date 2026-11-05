## Stereo perception: decide which chirality tags and double-bond cis/trans
## settings describe real stereo elements, and remove the rest.
##
## A tetrahedral centre is kept when its four positions (neighbours plus the
## implicit H or lone pair) are in distinct stereo-aware symmetry classes,
## or for ring stereo: two equivalent neighbours in a common ring that also
## holds another stereo candidate (cis/trans-1,4-disubstituted
## cyclohexanes). Double bonds keep cis/trans when they are outside rings
## smaller than 8 and each end has distinct substituents. Classes are
## recomputed after every removal, since dropping one element can make
## another non-stereogenic (or, through stereo labels, pseudoasymmetric
## centres distinguishable).

import std/[sets, sequtils]
import ./[mol, valence, rings, canon]

const
  minStereoRingSize = 8 ## double bonds in smaller rings are necessarily cis
  lonePairCentres = [15, 16, 33, 34]
    ## P, S, As, Se: three-coordinate centres whose lone pair holds the
    ## fourth position (phosphines, sulfoxides); N inverts so is excluded

func isBridgehead(i: int, rings, ringBondLists: seq[seq[int]]): bool =
  ## In two SSSR rings that share at least two bonds (norbornane-like).
  var mine: seq[int]
  for r, ring in rings:
    if i in ring:
      mine.add r
  for x in 0 ..< mine.len:
    for y in x + 1 ..< mine.len:
      if ringBondLists[mine[x]].countIt(it in ringBondLists[mine[y]]) >= 2:
        return true

func canBeTetrahedral(m: Mol, i: int, inThreeRing, bridgehead: bool): bool =
  ## Four stereo positions: four neighbours, three plus one H, or three
  ## plus a lone pair on P/S/As/Se. Three-coordinate N normally inverts,
  ## but not in a 3-membered ring or at a ring bridgehead (as in RDKit).
  let h = totalHCount(m, i)
  let d = m.degree(i)
  if d + h == 4 and h <= 1:
    return true
  if d != 3 or h != 0:
    return false
  let z = m.atoms[i].element
  if z in lonePairCentres:
    return true
  z == 7 and m.atoms[i].charge == 0 and (inThreeRing or bridgehead)

func canBeStereoDoubleBond(m: Mol, bi: int, smallRingBonds: HashSet[int]): bool =
  let bd = m.bonds[bi]
  if bd.order != boDouble or bi in smallRingBonds:
    return false
  for (endAtom, other) in [(bd.a, bd.b), (bd.b, bd.a)]:
    # one or two substituents; the other position may be an H or (on N) a
    # lone pair, as in oximes and imines
    let subs = substituents(m, endAtom, other)
    if subs.len notin 1 .. 2:
      return false
    if totalHCount(m, endAtom) + subs.len > 2:
      return false
    for s in subs:
      if m.bonds[m.bondBetween(endAtom, s)].order notin {boSingle, boAromatic}:
        return false
  true

func isPlainH(m: Mol, i: int): bool =
  ## An explicit hydrogen atom indistinguishable from an implicit one.
  let a = m.atoms[i]
  a.element == 1 and a.isotope == 0 and a.charge == 0 and a.atomMap == 0 and
    m.degree(i) == 1

func isRingStereo(
    m: Mol,
    i, twin: int,
    classes: seq[int],
    inRing: seq[bool],
    otherChiralInSystem: bool,
): bool =
  ## Ring stereo: the two equivalent neighbours (class of `twin`) are both
  ## ring neighbours, and the ring system holds another chiral atom.
  if not otherChiralInSystem:
    return false
  for (nb, bi) in m.adj[i]:
    if classes[nb] == classes[twin] and not inRing[bi]:
      return false
  true

proc assignStereo*(m: var Mol) =
  ## Removes chirality tags and double-bond stereo that do not describe a
  ## real stereo element (see the module docs).
  if not hasStereo(m):
    return # nothing to check; skip ring perception
  let rings = sssr(m)
  let ringBondLists = rings.mapIt(ringBondIndices(m, it))
  # bridgehead N: fixed by the ring system, so kept as a centre, but (in cage
  # aminals such as hexamethylenetetramine derivatives) not independent, so
  # excluded from the ring-stereo rule below
  var bridgeheadN = newSeq[bool](m.numAtoms)
  var anyStereo = false
  for i in 0 ..< m.numAtoms:
    if m.atoms[i].chiral == ctNone:
      continue
    let inThreeRing = rings.anyIt(it.len == 3 and i in it)
    bridgeheadN[i] = m.atoms[i].element == 7 and isBridgehead(i, rings, ringBondLists)
    if canBeTetrahedral(m, i, inThreeRing, bridgeheadN[i]):
      anyStereo = true
    else:
      m.atoms[i].chiral = ctNone
  var smallRingBonds = initHashSet[int]()
  for r, ring in rings:
    if ring.len < minStereoRingSize:
      for bi in ringBondLists[r]:
        smallRingBonds.incl bi
  for bi in 0 ..< m.numBonds:
    if m.bonds[bi].stereo != bsNone:
      if canBeStereoDoubleBond(m, bi, smallRingBonds):
        anyStereo = true
      else:
        m.clearDoubleBondStereo(bi)
  if not anyStereo:
    return

  # ring systems (components of the ring-bond graph), for the ring-stereo
  # rule; unlike a choice of SSSR rings these don't depend on atom order
  let inRing = ringBonds(m)
  let ringSystem = connectedComponents(m, inRing)

  # Each pass decides every removal from the state at the start of the pass
  # and applies them together, so the result doesn't depend on atom order.
  while true:
    let classes = symmetryClasses(m)
    var chiralInSystem = newSeq[int](m.numAtoms)
    for i in 0 ..< m.numAtoms:
      if m.atoms[i].chiral != ctNone:
        inc chiralInSystem[ringSystem[i]]
    var dropAtoms: seq[int]
    var dropBonds: seq[int]
    for i in 0 ..< m.numAtoms:
      if m.atoms[i].chiral == ctNone:
        continue
      var seen = initHashSet[int]()
      var dup: seq[int] # neighbours whose class was already seen
      for x in chiralReference(m, i):
        # an explicit [H] neighbour is the same as the implicit H slot
        let c =
          if x < 0 or isPlainH(m, x):
            -1
          else:
            classes[x]
        if c in seen:
          dup.add x
        seen.incl c
      if dup.len == 0:
        continue
      let keep =
        dup.len == 1 and dup[0] >= 0 and not isPlainH(m, dup[0]) and not bridgeheadN[i] and
        isRingStereo(m, i, dup[0], classes, inRing, chiralInSystem[ringSystem[i]] > 1)
      if not keep:
        dropAtoms.add i
    for bi in 0 ..< m.numBonds:
      let bd = m.bonds[bi]
      if bd.stereo == bsNone:
        continue
      for (endAtom, other) in [(bd.a, bd.b), (bd.b, bd.a)]:
        let subs = substituents(m, endAtom, other)
        if subs.len == 2 and classes[subs[0]] == classes[subs[1]]:
          dropBonds.add bi
          break
    if dropAtoms.len == 0 and dropBonds.len == 0:
      break
    for i in dropAtoms:
      m.atoms[i].chiral = ctNone
    for bi in dropBonds:
      m.clearDoubleBondStereo(bi)
