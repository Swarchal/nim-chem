## Explicit hydrogens: turning hydrogen atoms into H counts on their
## neighbours (`removeHs`) and back (`addHs`), keeping stereo.
##
## The SMILES and molfile readers run `removeHs` by default (as RDKit
## does), so `[H]C` and `C` give the same molecule. Hydrogen atoms that
## carry information are kept: isotopes (`[2H]`), charged or mapped
## hydrogens, H bonded to another H, a dummy atom or several atoms, and an
## H that is the only substituent defining a double bond's cis/trans
## (`[H]/N=C/C`) or that sits on a three-coordinate stereocentre whose lone
## pair holds the fourth position (`C[S@]([H])=O`).

import std/sequtils
import ./[mol, valence, geometry]

func isRemovableH*(m: Mol, i: int): bool =
  ## Whether hydrogen atom `i` can become an H count on its neighbour
  ## without losing information (see the module docs).
  let a = m.atoms[i]
  if a.element != 1 or a.isotope != 0 or a.charge != 0 or a.atomMap != 0 or
      m.degree(i) != 1:
    return false
  let (nb, bi) = m.adj[i][0]
  if m.atoms[nb].element <= 1 or m.bonds[bi].order != boSingle:
    return false
  # the only substituent fixing a double bond's cis/trans
  for (other, obi) in m.adj[nb]:
    let bd = m.bonds[obi]
    if bd.stereo != bsNone and i in [bd.stereoAtoms[0], bd.stereoAtoms[1]] and
        m.degree(nb) == 2:
      return false
  # a lone-pair stereocentre (three neighbours and no H) would lose a position
  if m.atoms[nb].chiral != ctNone and m.degree(nb) == 3 and totalHCount(m, nb) == 0:
    return false
  true

func removeHs*(m: Mol): Mol =
  ## A copy of `m` with the removable hydrogen atoms (`isRemovableH`)
  ## turned into H counts on their neighbours. Chirality tags and
  ## double-bond stereo are re-expressed; coordinates are kept for the
  ## remaining atoms.
  let remove = toSeq(0 ..< m.numAtoms).mapIt(isRemovableH(m, it))
  if not remove.anyIt(it):
    return m
  var src = m
  # double-bond references: switch from a removed H to the other substituent
  for bi, bd in m.bonds:
    if bd.stereo == bsNone:
      continue
    var refs = [bd.stereoAtoms[0], bd.stereoAtoms[1]]
    var cis = bd.stereo == bsCis
    for k, (endAtom, partner) in [(bd.a, bd.b), (bd.b, bd.a)]:
      if not remove[refs[k]]:
        continue
      let alt = toSeq(m.neighbors(endAtom)).filterIt(it notin [partner, refs[k]])
      if alt.len != 1 or remove[alt[0]]:
        refs[k] = -1
      else:
        refs[k] = alt[0]
        cis = not cis
    if refs[0] < 0 or refs[1] < 0:
      src.clearDoubleBondStereo(bi)
    else:
      src.setDoubleBondStereo(bi, refs[0], refs[1], cis)
  # H counts and chirality, worked out on the original numbering
  var hTarget = newSeq[int](m.numAtoms) # total H after removal, or -1
  var tags = newSeq[ChiralTag](m.numAtoms)
  for i in 0 ..< m.numAtoms:
    hTarget[i] = -1
    if remove[i]:
      continue
    var removed = 0
    var before: seq[int] # tag order with removed H as the H slot (-1)
    var after = @[-1] # the order once they are counts
    for x in chiralReference(m, i):
      if x >= 0 and remove[x]:
        inc removed
        before.add -1
      else:
        before.add x
        if x >= 0:
          after.add x
    if removed == 0:
      continue
    hTarget[i] = totalHCount(m, i) + removed
    tags[i] =
      if m.atoms[i].chiral == ctNone or before.countIt(it == -1) > 1:
        ctNone
      elif permutationParity(before, after) == 1:
        flip(m.atoms[i].chiral)
      else:
        m.atoms[i].chiral
  var keep: seq[int]
  for i in 0 ..< m.numAtoms:
    if not remove[i]:
      keep.add i
  result = subMol(src, keep)
  for ni, oi in keep:
    if hTarget[oi] < 0:
      continue
    if totalHCount(result, ni) != hTarget[oi]:
      result.atoms[ni].noImplicit = true
      result.atoms[ni].explicitH = hTarget[oi]
    result.atoms[ni].chiral = tags[oi]

func addHs*(m: Mol): Mol =
  ## A copy of `m` with every H count (explicit and implicit) turned into
  ## hydrogen atoms, appended after the existing atoms. Chirality tags are
  ## re-expressed; with coordinates, the new atoms are placed next to
  ## their neighbours (see `geometry.placeHydrogens`).
  result = m
  for i in 0 ..< m.numAtoms:
    let h = totalHCount(m, i)
    if h == 0:
      continue
    let before = chiralReference(m, i)
    result.atoms[i].explicitH = 0
    var hs: seq[int]
    for _ in 1 .. h:
      let hi = result.addAtom(initAtom(1))
      result.addBond(i, hi)
      hs.add hi
    if implicitHCount(result, i) != 0:
      result.atoms[i].noImplicit = true
    if m.atoms[i].chiral != ctNone:
      # the H slot (-1) becomes the first new H atom
      let mapped = before.mapIt(
        if it < 0:
          hs[0]
        else:
          it
      )
      let after = chiralReference(result, i)
      if mapped.len == after.len and permutationParity(mapped, after) == 1:
        result.atoms[i].chiral = flip(m.atoms[i].chiral)
    result.placeHydrogens(i, hs)
