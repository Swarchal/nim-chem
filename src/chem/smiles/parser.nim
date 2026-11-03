## SMILES parser (OpenSMILES subset).
##
## Supported: organic subset and bracket atoms (isotope, chirality `@`/`@@`
## and `@TH1`/`@TH2`, H count, charge, atom map), bonds `- = # $ : / \`,
## branches, ring closures (`0-9`, `%nn`) and `.` disconnections.

import std/[strutils, tables]
import ../[mol, elements, valence, rings]

type
  RingOpen = object
    atom: int
    order: BondOrder     ## boNone if unspecified
    stereo: BondStereo
    slot: int            ## position in nbrOrder[atom] reserved for the partner

  Parser = object
    s: string
    i: int
    m: Mol
    nbrOrder: seq[seq[int]]       ## neighbours in SMILES order; -1 = H
    openRings: Table[int, RingOpen]

proc fail(p: Parser; msg: string; pos = -1) {.noreturn.} =
  let at = if pos >= 0: pos else: p.i
  var e = newException(SmilesError, msg & " at position " & $at & " in '" & p.s & "'")
  e.pos = at
  raise e

func peek(p: Parser; off = 0): char =
  let j = p.i + off
  if j < p.s.len: p.s[j] else: '\0'

proc readInt(p: var Parser): int =
  ## Reads a run of decimal digits (returns -1 if none).
  if not p.peek.isDigit: return -1
  let start = p.i
  result = 0
  while p.peek.isDigit:
    if p.i - start >= 9: p.fail("number too large", start)
    result = result * 10 + (ord(p.peek) - ord('0'))
    inc p.i

proc addAtomNode(p: var Parser; a: Atom): int =
  result = p.m.addAtom(a)
  p.nbrOrder.add newSeq[int]()

proc connect(p: var Parser; a, b: int; order: BondOrder; stereo: BondStereo;
             pos: int) =
  ## Adds bond a->b. `order` == boNone means "unspecified"; it is resolved
  ## once the whole string has been read (see `resolveImplicitOrders`).
  try:
    p.m.addBond(a, b, order, stereo)
  except ChemError as e:
    p.fail(e.msg, pos)

proc parseOrganicAtom(p: var Parser): Atom =
  let c = p.peek
  var sym = $c
  if c == 'B' and p.peek(1) == 'r': sym = "Br"
  elif c == 'C' and p.peek(1) == 'l': sym = "Cl"
  if sym == "*":
    inc p.i
    return initAtom(0)
  if sym in organicSubset:
    p.i += sym.len
    return initAtom(atomicNumber(sym))
  if sym in organicAromatic:
    inc p.i
    return initAtom(atomicNumber(sym.toUpperAscii), aromatic = true)
  p.fail("unexpected character '" & c & "'")

proc parseBracketAtom(p: var Parser): Atom =
  let start = p.i
  inc p.i  # '['
  result.noImplicit = true
  let iso = p.readInt
  if iso > 0: result.isotope = iso
  # element symbol
  let c = p.peek
  if c == '*':
    result.element = 0
    inc p.i
  elif c.isLowerAscii:
    let two = c & p.peek(1)
    if two in aromaticSymbols:
      result.element = atomicNumber(two[0].toUpperAscii & two[1])
      p.i += 2
    elif $c in aromaticSymbols:
      result.element = atomicNumber($c.toUpperAscii)
      inc p.i
    else:
      p.fail("invalid aromatic symbol '" & c & "'")
    result.aromatic = true
  elif c.isUpperAscii:
    var z = -1
    if p.peek(1).isLowerAscii:
      z = atomicNumber(c & p.peek(1))
      if z >= 0: p.i += 2
    if z < 0:
      z = atomicNumber($c)
      if z < 0: p.fail("unknown element symbol")
      inc p.i
    result.element = z
  else:
    p.fail("expected element symbol in bracket atom")
  # chirality
  if p.peek == '@':
    inc p.i
    if p.peek == '@':
      inc p.i
      result.chiral = ctCW
    elif p.peek == 'T' and p.peek(1) == 'H':
      p.i += 2
      case p.readInt
      of 1: result.chiral = ctCCW
      of 2: result.chiral = ctCW
      else: p.fail("invalid @TH chirality class")
    elif p.peek.isUpperAscii and p.peek != 'H':   # @SP, @TB, @OH, @AL ...
      p.fail("only tetrahedral chirality (@, @@) is supported")
    else:
      result.chiral = ctCCW
  # hydrogens
  if p.peek == 'H':
    inc p.i
    let n = p.readInt
    result.explicitH = if n < 0: 1 else: n
  # charge
  if p.peek in {'+', '-'}:
    let ch = p.peek
    inc p.i
    var n = p.readInt
    if n < 0:   # "+", "++", "---" ...
      n = 1
      while p.peek == ch:
        inc n
        inc p.i
    result.charge = if ch == '+': n else: -n
  # atom map
  if p.peek == ':':
    inc p.i
    let n = p.readInt
    if n < 0: p.fail("expected atom map number after ':'")
    result.atomMap = n
  if p.peek != ']':
    p.fail("unterminated or malformed bracket atom", start)
  inc p.i

proc parseBondSymbol(c: char): (BondOrder, BondStereo) =
  case c
  of '-': (boSingle, bsNone)
  of '=': (boDouble, bsNone)
  of '#': (boTriple, bsNone)
  of '$': (boQuadruple, bsNone)
  of ':': (boAromatic, bsNone)
  of '/': (boSingle, bsUp)
  of '\\': (boSingle, bsDown)
  else: (boNone, bsNone)

proc fixChirality(p: var Parser) =
  ## Converts chirality from SMILES neighbour order to the Mol convention
  ## ([H] & adjacency order).
  for i in 0..<p.m.numAtoms:
    if p.m.atoms[i].chiral == ctNone: continue
    let reference = chiralReference(p.m, i)
    assert reference.len == p.nbrOrder[i].len
    if permutationParity(p.nbrOrder[i], reference) == 1:
      p.m.atoms[i].chiral = flip(p.m.atoms[i].chiral)

proc resolveImplicitOrders(m: var Mol) =
  ## An unwritten bond is aromatic if it joins two aromatic atoms within a
  ## ring, otherwise single (e.g. the biaryl bond in c1ccccc1c1ccccc1).
  let inRing = ringBonds(m)
  for bi, bd in m.bonds.mpairs:
    if bd.order == boNone:
      bd.order =
        if inRing[bi] and m.atoms[bd.a].aromatic and m.atoms[bd.b].aromatic: boAromatic
        else: boSingle

proc parseSmiles*(s: string; sanitize = true): Mol =
  ## Parses a SMILES string into a `Mol`. Raises `SmilesError` on syntax
  ## errors and (when `sanitize` is true) `ValenceError` on over-valent atoms.
  var p = Parser(s: s)
  var prev = -1
  var branchStack: seq[int]
  var pendingOrder = boNone
  var pendingStereo = bsNone
  var pendingPos = -1

  template clearPending() =
    pendingOrder = boNone
    pendingStereo = bsNone
    pendingPos = -1

  while p.i < s.len:
    let c = p.peek
    case c
    of '(':
      if prev < 0: p.fail("branch opened without a preceding atom")
      if pendingPos >= 0: p.fail("bond symbol before '('", pendingPos)
      branchStack.add prev
      inc p.i
    of ')':
      if branchStack.len == 0: p.fail("unmatched ')'")
      if pendingPos >= 0: p.fail("dangling bond symbol", pendingPos)
      prev = branchStack.pop
      inc p.i
    of '.':
      if pendingPos >= 0: p.fail("dangling bond symbol", pendingPos)
      prev = -1
      inc p.i
    of '-', '=', '#', '$', ':', '/', '\\':
      if pendingPos >= 0: p.fail("two consecutive bond symbols")
      (pendingOrder, pendingStereo) = parseBondSymbol(c)
      pendingPos = p.i
      inc p.i
    of '0'..'9', '%':
      if prev < 0: p.fail("ring closure without a preceding atom")
      let pos = p.i
      var num: int
      if c == '%':
        inc p.i
        if not (p.peek.isDigit and p.peek(1).isDigit):
          p.fail("'%' must be followed by two digits")
        num = (ord(p.peek) - ord('0')) * 10 + (ord(p.peek(1)) - ord('0'))
        p.i += 2
      else:
        num = ord(c) - ord('0')
        inc p.i
      var r: RingOpen
      if p.openRings.pop(num, r):
        if r.atom == prev: p.fail("ring closure to the same atom", pos)
        var order = r.order
        if pendingOrder != boNone:
          if order != boNone and order != pendingOrder:
            p.fail("conflicting ring closure bond orders", pos)
          order = pendingOrder
        # Each mark reads from its own atom, so consistent marks at the two
        # ends of a ring bond are opposite symbols (e.g. C/1 ... C\1).
        if r.stereo != bsNone and pendingStereo != bsNone and r.stereo == pendingStereo:
          p.fail("conflicting ring closure bond directions", pos)
        # store the bond in the direction its / or \ mark was written
        if r.stereo != bsNone or pendingStereo == bsNone:
          p.connect(r.atom, prev, order, r.stereo, pos)
        else:
          p.connect(prev, r.atom, order, pendingStereo, pos)
        p.nbrOrder[r.atom][r.slot] = prev
        p.nbrOrder[prev].add r.atom
      else:
        p.openRings[num] = RingOpen(atom: prev, order: pendingOrder,
                                    stereo: pendingStereo,
                                    slot: p.nbrOrder[prev].len)
        p.nbrOrder[prev].add -2   # placeholder until the ring closes
      clearPending()
    else:
      let a =
        if c == '[': p.parseBracketAtom()
        else: p.parseOrganicAtom()
      let idx = p.addAtomNode(a)
      if prev >= 0:
        p.connect(prev, idx, pendingOrder, pendingStereo, pendingPos)
        p.nbrOrder[prev].add idx
        p.nbrOrder[idx].add prev
      elif pendingPos >= 0:
        p.fail("bond symbol without a preceding atom", pendingPos)
      if a.explicitH > 0:
        p.nbrOrder[idx].add -1
      clearPending()
      prev = idx

  if pendingPos >= 0: p.fail("dangling bond symbol at end of input", pendingPos)
  if branchStack.len > 0: p.fail("unclosed branch '('")
  for num in p.openRings.keys:
    p.fail("unclosed ring bond " & $num)

  resolveImplicitOrders(p.m)

  p.fixChirality()
  if sanitize: valence.sanitize(p.m)
  p.m
