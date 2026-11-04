## Valence model: implicit hydrogen counts and valence validation.
##
## Allowed valences come from `elements.defaultValences`, adjusted for charge
## using the isoelectronic element (e.g. N+ behaves like C, O- like F).
## Aromatic bonds count as 1; an aromatic atom gains one extra valence unit
## (from the pi system) when that still fits its lowest default valence. This
## approximates the Kekulé form without performing Kekulization.

import ./[mol, elements]

func allowedValences*(a: Atom): seq[int] =
  ## Allowed total valences for an atom, accounting for formal charge.
  ## Empty if the element has no default valence model (e.g. most metals).
  if elementTable[a.element].valences.len == 0: return @[]
  let z = a.element - a.charge
  if z in 1..maxAtomicNumber: elementTable[z].valences else: @[]

func explicitValence*(m: Mol; i: int): int =
  ## Sum of bond valences to neighbouring atoms (aromatic bonds count 1).
  for (_, bi) in m.adj[i]:
    result += bondValence(m.bonds[bi].order)

func hasAromaticBond(m: Mol; i: int): bool =
  for (_, bi) in m.adj[i]:
    if m.bonds[bi].order == boAromatic: return true

func defaultImplicitHCount*(m: Mol; i: int): int =
  ## Number of hydrogens a SMILES reader would assign to atom `i` if it were
  ## written as an organic-subset (non-bracket) atom.
  let a = m.atoms[i]
  let vals = allowedValences(a)
  if vals.len == 0: return 0
  let v = explicitValence(m, i)
  if a.aromatic and hasAromaticBond(m, i):
    let t = vals[0]
    if v + 1 <= t: return t - (v + 1)
    if v <= t: return t - v
    for val in vals:
      if val >= v + 1: return val - (v + 1)
    return 0
  for val in vals:
    if val >= v: return val - v
  0

func implicitHCount*(m: Mol; i: int): int =
  ## Implicit hydrogens on atom `i` (0 for bracket atoms / `noImplicit`).
  if m.atoms[i].noImplicit: 0 else: defaultImplicitHCount(m, i)

func totalHCount*(m: Mol; i: int; includeNeighbors = false): int =
  ## Hydrogens attached to atom `i`: explicit count + implicit. With
  ## `includeNeighbors`, also counts bonded hydrogen atoms.
  result = m.atoms[i].explicitH + implicitHCount(m, i)
  if includeNeighbors:
    for nb in m.neighbors(i):
      if m.atoms[nb].element == 1: inc result

func totalValence*(m: Mol; i: int): int =
  explicitValence(m, i) + m.atoms[i].explicitH + implicitHCount(m, i)

func validateValences*(m: Mol): seq[string] =
  ## Returns a description of every over-valent atom (empty if all are fine).
  for i, a in m.atoms:
    let vals = allowedValences(a)
    if vals.len == 0: continue
    let total = totalValence(m, i)
    if total > vals[^1]:
      result.add "atom " & $i & " (" & $a & ") has valence " & $total &
                 ", maximum allowed is " & $vals[^1]

proc checkValences*(m: Mol) =
  ## Raises `ValenceError` if any atom is over-valent.
  let errs = validateValences(m)
  if errs.len > 0:
    raise newException(ValenceError, "invalid valence: " & errs[0])

func chiralReference*(m: Mol; i: int): seq[int] =
  ## Neighbour order that atom `i`'s chirality tag refers to: -1 (the
  ## hydrogen count) first if the atom has any H, then `adj` order. For a
  ## three-coordinate centre (sulfoxide, phosphine) the lone pair has a
  ## fixed position and is not listed, matching RDKit: C[S@](=O)c1ccccc1 is
  ## the same as [S@](C)(=O)c1ccccc1.
  if totalHCount(m, i) > 0: result.add -1
  for nb in m.neighbors(i): result.add nb
