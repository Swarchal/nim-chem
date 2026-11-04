## Stereo from atom coordinates, and back: tetrahedral chirality from 3D
## positions or 2D wedge/hash bonds, double-bond cis/trans from either, and
## MDL atom parities (which depend only on atom numbering).
##
## A chirality tag (see `mol`) describes the neighbours in
## `chiralReference` order: looking from the first towards the centre, the
## others run anticlockwise (`ctCCW`, SMILES `@`) or clockwise. An implicit
## H comes first; a lone pair (three-coordinate S, P, N) is not listed, and
## counts as last (see `stereoSlots`), as RDKit places it.

import std/[algorithm, math, sequtils]
import ./[mol, valence, rings]

type BondWedge* = enum
  ## a 2D stereo bond, whose narrow end is at `bond.a`
  bwNone, bwWedge, bwHash

func `-`(a, b: Point3): Point3 = (a.x - b.x, a.y - b.y, a.z - b.z)
func `+`(a, b: Point3): Point3 = (a.x + b.x, a.y + b.y, a.z + b.z)
func `*`(a: Point3; k: float): Point3 = (a.x * k, a.y * k, a.z * k)
func dot(a, b: Point3): float = a.x * b.x + a.y * b.y + a.z * b.z
func cross(a, b: Point3): Point3 =
  (a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x)
func norm(a: Point3): float = sqrt(dot(a, a))

func is3D*(m: Mol): bool =
  ## Whether `m` has coordinates that are not all in the z = 0 plane.
  for p in m.coords:
    if p.z != 0.0: return true

func meanBondLength(m: Mol): float =
  for bd in m.bonds: result += norm(m.coords[bd.a] - m.coords[bd.b])
  if m.numBonds > 0: result /= m.numBonds.float

func stereoSlots(m: Mol; i: int): seq[int] =
  ## The four stereo positions of atom `i` in the order its tag refers to:
  ## `chiralReference` (an H slot, -1, first), then for three neighbours
  ## with no H a lone-pair slot, -2, last (where RDKit puts it). Empty if
  ## the atom doesn't have four positions.
  let h = totalHCount(m, i)
  let d = m.degree(i)
  if d + h == 4 and h <= 1: chiralReference(m, i)
  elif d == 3 and h == 0: chiralReference(m, i) & @[-2]
  else: @[]

func chiralFromPoints(c: Point3; pts: openArray[Point3]): ChiralTag =
  ## The tag for 4 neighbour positions around centre `c`, or `ctNone` if
  ## they are (nearly) coplanar: the signed volume of the neighbours is
  ## positive when, seen from the first, the others run anticlockwise.
  let t = dot(cross(pts[2] - pts[1], pts[3] - pts[1]), pts[0] - pts[1])
  var scale = 0.0
  for p in pts: scale = max(scale, norm(p - c))
  if abs(t) <= 1e-3 * scale * scale * scale: ctNone
  elif t > 0: ctCCW
  else: ctCW

func tagFromPositions(m: Mol; i: int; slots: openArray[int];
                      lift: openArray[float]): ChiralTag =
  ## Atom `i`'s tag from the coordinates, with neighbour `nb` raised by
  ## `lift[nb]` in z (for wedges); an implicit H or lone pair is placed
  ## opposite the neighbours.
  template pos(nb: int): Point3 =
    (m.coords[nb].x, m.coords[nb].y, m.coords[nb].z + lift[nb])
  let c = m.coords[i]
  var sum: Point3 = (0.0, 0.0, 0.0)
  for nb in m.neighbors(i):
    let d = pos(nb) - c
    let l = norm(d)
    if l == 0: return ctNone
    sum = sum + d * (1 / l)
  var pts: array[4, Point3]
  for k, x in slots: pts[k] = if x < 0: c - sum else: pos(x)
  chiralFromPoints(c, pts)

proc stereoFromCoords*(m: var Mol; wedges: openArray[BondWedge];
                       unknownBonds: openArray[bool]) =
  ## Sets chirality tags and double-bond cis/trans from `m.coords`. In 3D
  ## every atom with four stereo positions gets a tag; in 2D only atoms at
  ## the narrow end of a wedge or hash bond (and none with a wavy bond,
  ## `unknownBonds`, there). Double bonds marked in `unknownBonds`
  ## ("either"), or next to a wavy single bond, get no stereo. Run
  ## `stereo.assignStereo` (part of `sanitize`) afterwards to remove tags
  ## that don't describe real stereo.
  if m.coords.len == 0: return
  let threeD = m.is3D
  let bondLen = meanBondLength(m)
  var lift = newSeq[float](m.numAtoms)   # z offsets from wedges at one atom
  for i in 0..<m.numAtoms:
    m.atoms[i].chiral = ctNone
    let slots = stereoSlots(m, i)
    if slots.len == 0: continue
    var wedged, wavy = false
    for (nb, bi) in m.adj[i]:
      lift[nb] = 0
      if m.bonds[bi].a != i: continue
      if unknownBonds[bi]: wavy = true
      if wedges[bi] == bwNone: continue
      lift[nb] = if wedges[bi] == bwWedge: bondLen else: -bondLen
      wedged = true
    if wavy or not threeD and not wedged: continue
    m.atoms[i].chiral = tagFromPositions(m, i, slots, lift)
  for bi in 0..<m.numBonds:
    m.clearDoubleBondStereo(bi)
    let bd = m.bonds[bi]
    if bd.order != boDouble or unknownBonds[bi]: continue
    var refs: array[2, int]
    var ok = true
    for k, (endAtom, other) in [(bd.a, bd.b), (bd.b, bd.a)]:
      refs[k] = -1
      for (nb, nbi) in m.adj[endAtom]:
        if nb == other: continue
        if unknownBonds[nbi]: ok = false
        if refs[k] < 0: refs[k] = nb
      if refs[k] < 0: ok = false
    if not ok: continue
    let axis = m.coords[bd.b] - m.coords[bd.a]
    let l = norm(axis)
    if l == 0: continue
    let u = axis * (1 / l)
    template perp(v: Point3): Point3 = v - u * dot(v, u)
    let pa = perp(m.coords[refs[0]] - m.coords[bd.a])
    let pb = perp(m.coords[refs[1]] - m.coords[bd.b])
    let d = dot(pa, pb)
    if abs(d) <= 1e-3 * norm(pa) * norm(pb):
      continue   # a substituent in line with the bond
    m.setDoubleBondStereo(bi, refs[0], refs[1], cis = d > 0)

# ---------------------------------------------------------------------------
# MDL parity

func chiralFromParity*(m: Mol; i, parity: int): ChiralTag =
  ## Atom `i`'s tag from an MDL parity (1 or 2; anything else: `ctNone`).
  ## Parity 1: with the neighbours numbered in atom order (an implicit H or
  ## lone pair last), 1-2-3 run clockwise seen with 4 pointing away, i.e.
  ## anticlockwise seen from 4: the tag `ctCCW` for (4, 1, 2, 3).
  let slots = stereoSlots(m, i)
  if parity notin 1..2 or slots.len == 0: return ctNone
  var numbered = slots
  numbered.sort(proc (a, b: int): int =
    if a < 0: 1 elif b < 0: -1 else: cmp(a, b))
  let viewFrom4 = @[numbered[3], numbered[0], numbered[1], numbered[2]]
  if (permutationParity(viewFrom4, slots) == 1) xor (parity == 2): ctCW
  else: ctCCW

func mdlParity*(m: Mol; i: int): int =
  ## The MDL parity (1 or 2) of atom `i`'s chirality tag, or 0.
  let tag = m.atoms[i].chiral
  if tag == ctNone or stereoSlots(m, i).len == 0: 0
  elif chiralFromParity(m, i, 1) == tag: 1
  else: 2

# ---------------------------------------------------------------------------
# Choosing wedges

proc chooseWedges*(m: Mol): seq[(int, BondWedge)] =
  ## For 2D coordinates: per bond, the atom at the narrow end and a wedge
  ## or hash that gives each chiral atom its tag (`(-1, bwNone)` for plain
  ## bonds, and for every bond without 2D coordinates). Prefers bonds to
  ## non-chiral, non-ring, low-degree neighbours; the atoms with the fewest
  ## usable bonds choose first. Raises `ChemError` if a chiral atom can't
  ## get a wedge (every single bond at it is needed by another chiral atom,
  ## or the drawing can't show its configuration).
  result = newSeqWith(m.numBonds, (-1, bwNone))
  if m.coords.len == 0 or m.is3D or m.atoms.allIt(it.chiral == ctNone): return
  let inRing = ringBonds(m)
  let bondLen = meanBondLength(m)
  var lift = newSeq[float](m.numAtoms)
  # per chiral atom, the bonds that could show its tag, best first
  var options: seq[(int, seq[(int, BondWedge)])]   # (atom, [(bond, kind)])
  for i in 0..<m.numAtoms:
    let tag = m.atoms[i].chiral
    let slots = stereoSlots(m, i)
    if tag == ctNone or slots.len == 0: continue
    var ranked: seq[(int, int, int, int)]   # (chiral nb, ring, degree, bond)
    for (nb, bi) in m.adj[i]:
      if m.bonds[bi].order == boSingle:
        ranked.add (ord(m.atoms[nb].chiral != ctNone), ord(inRing[bi]), m.degree(nb), bi)
    ranked.sort()
    var usable: seq[(int, BondWedge)]
    for (_, _, _, bi) in ranked:
      let up = m.bonds[bi].otherAtom(i)
      lift[up] = bondLen
      let t = tagFromPositions(m, i, slots, lift)
      lift[up] = 0
      if t != ctNone: usable.add (bi, if t == tag: bwWedge else: bwHash)
    options.add (i, usable)
  options.sort(proc (x, y: (int, seq[(int, BondWedge)])): int =
    cmp(x[1].len, y[1].len))
  var owner = newSeqWith(m.numBonds, -1)   # option index using each bond
  template free(k: int): int =
    ## the first option of atom `k` whose bond is free, or -1
    var found = -1
    for o, opt in options[k][1]:
      if owner[opt[0]] < 0:
        found = o
        break
    found
  template take(k, o: int) =
    let (bi, kind) = options[k][1][o]
    owner[bi] = k
    result[bi] = (options[k][0], kind)
  for k in 0..<options.len:
    var o = free(k)
    if o < 0:
      # take a bond from an atom that has another free one
      for o2, (bi, _) in options[k][1]:
        let j = owner[bi]
        let alt = free(j)
        if alt >= 0:
          result[bi] = (-1, bwNone)
          owner[bi] = -1
          take(j, alt)
          o = o2
          break
    if o < 0:
      raise newException(ChemError, "can't draw a wedge for chiral atom " &
                         $options[k][0] & " in these 2D coordinates")
    take(k, o)
