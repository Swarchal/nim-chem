## Core molecule graph types.
##
## A `Mol` is a value type holding atoms and bonds in flat sequences, with an
## adjacency list of `(neighbor, bondIndex)` pairs per atom. Hydrogens may be
## stored either as explicit atoms or as counts on heavy atoms (`explicitH`
## plus implicit hydrogens computed by the valence model).
##
## Chirality convention: an atom's `chiral` tag describes the order
##   [implicit slot (if any)] & neighbors in `adj[i]` order,
## looking from the first neighbour: `ctCCW` = remaining neighbours are
## anticlockwise (SMILES `@`), `ctCW` = clockwise (`@@`). The implicit slot
## is the atom's hydrogen count, if it has any (see
## `valence.chiralReference`).
##
## Double-bond stereo is stored on the double bond as cis/trans relative to
## one reference neighbour of each end (`stereoAtoms`); SMILES `/` `\` marks
## are only a syntax used by the parser and writer.

import std/[strutils, algorithm]
import ./elements

type
  ChemError* = object of CatchableError
  SmilesError* = object of ChemError
    pos*: int   ## 0-based character offset where the error was detected
  ValenceError* = object of ChemError

  BondOrder* = enum
    boNone, boSingle, boDouble, boTriple, boQuadruple, boAromatic

  BondStereo* = enum
    bsNone   ## unspecified
    bsCis    ## the two `stereoAtoms` are on the same side (Z-like)
    bsTrans  ## the two `stereoAtoms` are on opposite sides (E-like)

  ChiralTag* = enum
    ctNone, ctCCW, ctCW

  Atom* = object
    element*: int       ## atomic number (0 = wildcard '*')
    charge*: int
    isotope*: int       ## mass number, 0 = natural abundance
    explicitH*: int     ## hydrogens stored as a count on this atom
    noImplicit*: bool   ## if true, no implicit hydrogens are added (bracket atoms)
    aromatic*: bool
    chiral*: ChiralTag
    atomMap*: int

  Bond* = object
    a*, b*: int
    order*: BondOrder
    stereo*: BondStereo      ## double bonds only
    stereoAtoms*: (int, int) ## reference neighbours of `a` and `b` for `stereo`

  Mol* = object
    atoms*: seq[Atom]
    bonds*: seq[Bond]
    adj*: seq[seq[(int, int)]]  ## per atom: (neighbor atom, bond index)

func initAtom*(element: int; charge = 0; isotope = 0; explicitH = 0;
               noImplicit = false; aromatic = false; chiral = ctNone;
               atomMap = 0): Atom =
  Atom(element: element, charge: charge, isotope: isotope, explicitH: explicitH,
       noImplicit: noImplicit, aromatic: aromatic, chiral: chiral, atomMap: atomMap)

func numAtoms*(m: Mol): int {.inline.} = m.atoms.len
func numBonds*(m: Mol): int {.inline.} = m.bonds.len
func degree*(m: Mol; i: int): int {.inline.} = m.adj[i].len

func addAtom*(m: var Mol; a: Atom): int {.discardable.} =
  ## Appends an atom and returns its index.
  if not isValidElement(a.element):
    raise newException(ChemError, "invalid atomic number: " & $a.element)
  m.atoms.add a
  m.adj.add newSeq[(int, int)]()
  m.atoms.high

func bondBetween*(m: Mol; a, b: int): int =
  ## Index of the bond joining atoms `a` and `b`, or -1.
  for (nb, bi) in m.adj[a]:
    if nb == b: return bi
  -1

func addBond*(m: var Mol; a, b: int; order = boSingle): int {.discardable.} =
  ## Adds a bond between atoms `a` and `b`, returning its index.
  if a notin 0..m.atoms.high or b notin 0..m.atoms.high:
    raise newException(ChemError, "bond references a non-existent atom")
  if a == b:
    raise newException(ChemError, "cannot bond atom " & $a & " to itself")
  if m.bondBetween(a, b) >= 0:
    raise newException(ChemError, "atoms " & $a & " and " & $b & " are already bonded")
  m.bonds.add Bond(a: a, b: b, order: order)
  let bi = m.bonds.high
  m.adj[a].add (b, bi)
  m.adj[b].add (a, bi)
  bi

func otherAtom*(bd: Bond; i: int): int {.inline.} =
  if bd.a == i: bd.b else: bd.a

func setDoubleBondStereo*(m: var Mol; bi, refA, refB: int; cis: bool) =
  ## Sets cis/trans stereo on double bond `bi`: `refA` must be a neighbour
  ## of `bonds[bi].a` and `refB` of `bonds[bi].b` (other than each other).
  let bd = m.bonds[bi]
  if bd.order != boDouble:
    raise newException(ChemError, "bond " & $bi & " is not a double bond")
  if refA == bd.b or m.bondBetween(bd.a, refA) < 0 or
     refB == bd.a or m.bondBetween(bd.b, refB) < 0:
    raise newException(ChemError, "stereo reference atoms must be neighbours of the bond ends")
  m.bonds[bi].stereo = if cis: bsCis else: bsTrans
  m.bonds[bi].stereoAtoms = (refA, refB)

func clearDoubleBondStereo*(m: var Mol; bi: int) =
  m.bonds[bi].stereo = bsNone
  m.bonds[bi].stereoAtoms = (0, 0)

iterator neighbors*(m: Mol; i: int): int =
  for (nb, _) in m.adj[i]: yield nb

iterator atomBonds*(m: Mol; i: int): (int, int) =
  ## Yields `(neighbor, bondIndex)` pairs for atom `i`.
  for p in m.adj[i]: yield p

func bondValence*(o: BondOrder): int =
  ## Integer contribution of a bond to valence (aromatic counts as 1).
  case o
  of boNone: 0
  of boSingle, boAromatic: 1
  of boDouble: 2
  of boTriple: 3
  of boQuadruple: 4

func permutationParity*(a, b: openArray[int]): int =
  ## Parity (0 even, 1 odd) of the permutation taking sequence `a` to `b`.
  ## Both must contain the same distinct elements.
  assert a.len == b.len
  var p = newSeq[int](a.len)
  for i, x in a:
    let j = b.find(x)
    assert j >= 0, "sequences differ"
    p[i] = j
  var inv = 0
  for i in 0..<p.len:
    for j in i+1..<p.len:
      if p[i] > p[j]: inc inv
  inv and 1

func flip*(c: ChiralTag): ChiralTag =
  case c
  of ctNone: ctNone
  of ctCCW: ctCW
  of ctCW: ctCCW

func subMol*(m: Mol; atoms: openArray[int]): Mol =
  ## The molecule formed by `atoms` (in that order) and the bonds between
  ## them, which keep their relative order. Chirality tags are re-expressed
  ## for the new neighbour order; stereo is dropped where an atom it refers
  ## to is left out. Raises `ChemError` if an index is out of range or
  ## repeated.
  var oldToNew = newSeq[int](m.atoms.len)
  oldToNew.fill(-1)
  for ni, oi in atoms:
    if oi notin 0..m.atoms.high:
      raise newException(ChemError, "atom index " & $oi & " out of range")
    if oldToNew[oi] >= 0:
      raise newException(ChemError, "atom " & $oi & " listed twice")
    oldToNew[oi] = ni
    result.addAtom m.atoms[oi]
  for bd in m.bonds:
    let (na, nb) = (oldToNew[bd.a], oldToNew[bd.b])
    if na < 0 or nb < 0: continue
    let bi = result.addBond(na, nb, bd.order)
    if bd.stereo != bsNone and oldToNew[bd.stereoAtoms[0]] >= 0 and
       oldToNew[bd.stereoAtoms[1]] >= 0:
      result.bonds[bi].stereo = bd.stereo
      result.bonds[bi].stereoAtoms = (oldToNew[bd.stereoAtoms[0]], oldToNew[bd.stereoAtoms[1]])
  for ni, oi in atoms:
    if m.atoms[oi].chiral == ctNone: continue
    var before: seq[int]
    for (nb, _) in m.adj[oi]: before.add oldToNew[nb]
    var after: seq[int]
    for (nb, _) in result.adj[ni]: after.add nb
    if before.len != after.len:
      result.atoms[ni].chiral = ctNone   # a neighbour was left out
    elif permutationParity(before, after) == 1:
      result.atoms[ni].chiral = flip(result.atoms[ni].chiral)

func renumberAtoms*(m: Mol; newOrder: openArray[int]): Mol =
  ## Returns a copy with atoms reordered: new atom `i` is old atom `newOrder[i]`.
  ## Bond indices are unchanged and stereochemistry is preserved. The H slot
  ## of the chirality reference order (see module docs) is first both before
  ## and after renumbering, so only the neighbour order matters.
  ## Raises `ChemError` if `newOrder` is not a permutation.
  if newOrder.len != m.atoms.len:
    raise newException(ChemError, "newOrder must have one entry per atom")
  subMol(m, newOrder)   # rejects repeated and out-of-range indices

func chargeSuffix*(q: int): string =
  ## SMILES-style charge text: "", "+", "-", "+2", "-3", ...
  if q == 0: return ""
  result = if q > 0: "+" else: "-"
  if abs(q) > 1: result &= $abs(q)

func `$`*(a: Atom): string =
  result = symbol(a.element)
  if a.aromatic: result = result.toLowerAscii
  if a.isotope != 0: result = $a.isotope & result
  result &= chargeSuffix(a.charge)
  if a.explicitH > 0: result &= " H" & $a.explicitH
  if a.chiral != ctNone: result &= " " & $a.chiral
  if a.atomMap != 0: result &= " :" & $a.atomMap

func `$`*(m: Mol): string =
  result = "Mol(" & $m.numAtoms & " atoms, " & $m.numBonds & " bonds)\n"
  for i, a in m.atoms:
    result &= "  atom " & $i & ": " & $a & "\n"
  for i, bd in m.bonds:
    result &= "  bond " & $i & ": " & $bd.a & "-" & $bd.b & " " & $bd.order
    if bd.stereo != bsNone:
      result &= " " & $bd.stereo & " (" & $bd.stereoAtoms[0] & ", " & $bd.stereoAtoms[1] & ")"
    result &= "\n"
