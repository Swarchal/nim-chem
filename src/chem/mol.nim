## Core molecule graph types.
##
## A `Mol` is a value type holding atoms and bonds in flat sequences, with an
## adjacency list of `(neighbor, bondIndex)` pairs per atom. Hydrogens may be
## stored either as explicit atoms or as counts on heavy atoms (`explicitH`
## plus implicit hydrogens computed by the valence model).
##
## Chirality convention: an atom's `chiral` tag describes the order
##   [implicit/explicit-count H (if any)] & neighbors in `adj[i]` order,
## looking from the first neighbour: `ctCCW` = remaining neighbours are
## anticlockwise (SMILES `@`), `ctCW` = clockwise (`@@`).

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
    bsNone   ## no directional marker
    bsUp     ## '/' when read from atom `a` to atom `b`
    bsDown   ## '\' when read from atom `a` to atom `b`

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
    stereo*: BondStereo

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

func addBond*(m: var Mol; a, b: int; order = boSingle;
              stereo = bsNone): int {.discardable.} =
  ## Adds a bond between atoms `a` and `b`, returning its index.
  if a notin 0..m.atoms.high or b notin 0..m.atoms.high:
    raise newException(ChemError, "bond references a non-existent atom")
  if a == b:
    raise newException(ChemError, "cannot bond atom " & $a & " to itself")
  if m.bondBetween(a, b) >= 0:
    raise newException(ChemError, "atoms " & $a & " and " & $b & " are already bonded")
  m.bonds.add Bond(a: a, b: b, order: order, stereo: stereo)
  let bi = m.bonds.high
  m.adj[a].add (b, bi)
  m.adj[b].add (a, bi)
  bi

func otherAtom*(bd: Bond; i: int): int {.inline.} =
  if bd.a == i: bd.b else: bd.a

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

func renumberAtoms*(m: Mol; newOrder: openArray[int]): Mol =
  ## Returns a copy with atoms reordered: new atom `i` is old atom `newOrder[i]`.
  ## Chirality tags are adjusted so that the stereochemistry is preserved.
  ## The H slot of the chirality reference order (see module docs) is first
  ## both before and after renumbering, so only the neighbour order matters.
  if newOrder.len != m.atoms.len:
    raise newException(ChemError, "newOrder must have one entry per atom")
  var oldToNew = newSeq[int](m.atoms.len)
  oldToNew.fill(-1)
  for ni, oi in newOrder:
    if oi notin 0..m.atoms.high or oldToNew[oi] >= 0:
      raise newException(ChemError, "newOrder is not a permutation")
    oldToNew[oi] = ni
  for oi in newOrder:
    result.addAtom m.atoms[oi]
  # keep bonds in original order so bond indices are stable
  for bd in m.bonds:
    result.addBond(oldToNew[bd.a], oldToNew[bd.b], bd.order, bd.stereo)
  for ni, oi in newOrder:
    if m.atoms[oi].chiral != ctNone:
      var before: seq[int]
      for (nb, _) in m.adj[oi]: before.add oldToNew[nb]
      var after: seq[int]
      for (nb, _) in result.adj[ni]: after.add nb
      if permutationParity(before, after) == 1:
        result.atoms[ni].chiral = flip(result.atoms[ni].chiral)

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
    if bd.stereo != bsNone: result &= " " & $bd.stereo
    result &= "\n"
