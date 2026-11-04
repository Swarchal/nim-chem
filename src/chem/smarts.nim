## SMARTS parser (Daylight SMARTS subset).
##
## Atoms: organic-subset symbols, `*`, `a`, `A` outside brackets; inside
## brackets the primitives `*` `a` `A` `#n` element symbols (uppercase
## aliphatic, lowercase aromatic) `Dn` `Hn` `hn` `Rn` `rn` `vn` `Xn` `xn`
## `+n` `-n` isotope `@` `@@` `:n` (atom map) and recursive `$(...)`,
## combined with `!` `&` `,` `;` (highest to lowest precedence; adjacent
## primitives are joined by `&`).
##
## Bonds: `-` `=` `#` `$` `:` `~` `@` `/` `\` with the same operators; an
## unwritten bond means single or aromatic. Branches, ring closures and `.`
## work as in SMILES. Component-level grouping (`(C).(C)`) and reaction
## SMARTS are not supported.

import std/[strutils, tables, sequtils]
import ./[mol, elements, query]
import ./smiles/directions

type
  RingOpen = object
    atom: int
    bond: Expr[BondPrim]   ## nil if unwritten
    bondText: string       ## the bond expression as written
    dir: Dir
    slot: int              ## position in nbrOrder[atom] reserved for the partner

  Parser = object
    s: string
    i: int
    offset: int            ## position of `s` in the top-level string (for errors)
    q: Query
    marks: seq[Dir]        ## per bond: a `/` `\` mark as written from `bond.a`
    bondPos: seq[int]      ## per bond: input position (for errors)
    openRings: Table[int, RingOpen]
    nRecursive: int        ## recursive query ids handed out so far (shared
                           ## with sub-parsers)
    atom: AtomState        ## the bracket atom being read

  AtomState = object
    ## per bracket atom: what primitives set besides the expression
    chiral: ChiralTag
    hasH: bool       ## an H count > 0 (gives the chirality an H slot)
    atomMap: int
    first: bool      ## no primitive read yet (for `[H]` vs `[CH]`)

proc fail(p: Parser; msg: string; pos = -1) {.noreturn.} =
  let at = (if pos >= 0: pos else: p.i) + p.offset
  var e = newException(SmartsError, msg & " at position " & $at & " in SMARTS")
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

proc readIntOr(p: var Parser; default: int): int =
  result = p.readInt
  if result < 0: result = default

func elementExpr(z: int; aromatic: bool): Expr[AtomPrim] =
  andE(atomPrim(apAtomicNum, z), atomPrim(if aromatic: apAromatic else: apAliphatic))

proc parseSmartsInto(p: var Parser)
proc assignDoubleBondStereo(p: var Parser)

# ---------------------------------------------------------------------------
# Expressions

proc parseExpr[P](p: var Parser; parsePrim: proc (p: var Parser): Expr[P] {.nimcall.};
                  continues: proc (c: char): bool {.nimcall.}; always: P): Expr[P] =
  ## `;` (lowest precedence) over `,` over `&` or juxtaposition over `!`.
  ## `parsePrim` may return nil for tokens that only set parser state;
  ## `continues(c)` says whether `c` can follow within the same `&` term.
  var low: seq[Expr[P]]
  while true:
    var ors: seq[Expr[P]]
    while true:
      var ands: seq[Expr[P]]
      while true:
        var negate = false
        while p.peek == '!':
          negate = not negate
          inc p.i
        let e = parsePrim(p)
        if not e.isNil: ands.add(if negate: notE(e) else: e)
        if p.peek == '&': inc p.i
        elif not continues(p.peek): break
      ors.add(if ands.len == 0: prim(always) else: andE(ands))
      if p.peek != ',': break
      inc p.i
    low.add orE(ors)
    if p.peek != ';': break
    inc p.i
  andE(low)

# ---------------------------------------------------------------------------
# Atom expressions

proc parseRecursive(p: var Parser): Expr[AtomPrim] =
  ## `$(` ... `)`; `p.i` is at `$`.
  let start = p.i
  p.i += 2
  var depth = 1
  var j = p.i
  while j < p.s.len and depth > 0:
    if p.s[j] == '(': inc depth
    elif p.s[j] == ')': dec depth
    if depth > 0: inc j
  if depth > 0: p.fail("unterminated recursive SMARTS '$('", start)
  var sub = Parser(s: p.s[p.i ..< j], offset: p.offset + p.i,
                   nRecursive: p.nRecursive)
  if sub.s.len == 0: p.fail("empty recursive SMARTS", start)
  sub.parseSmartsInto()
  p.nRecursive = sub.nRecursive
  p.i = j + 1
  var r = new Query
  r[] = sub.q
  result = prim(AtomPrim(kind: apRecursive, sub: r, subId: p.nRecursive))
  inc p.nRecursive

proc countPrim(p: var Parser; kind: AtomPrimKind): Expr[AtomPrim] =
  ## `h`, `x`: a count, or with no number at least one.
  let n = p.readInt
  if n < 0: notE(atomPrim(kind, 0)) else: atomPrim(kind, n)

proc ringPrim(p: var Parser; kind: AtomPrimKind): Expr[AtomPrim] =
  ## `R`, `r`: with no number in a ring, `0` not in a ring, else `kind`.
  let n = p.readInt
  if n < 0: atomPrim(apInRing)
  elif n == 0: notE(atomPrim(apInRing))
  else: atomPrim(kind, n)

proc parseAtomPrim(p: var Parser): Expr[AtomPrim] =
  ## One atom primitive; nil for `@`, `@@` and `:n`, which only set state.
  let c = p.peek
  let start = p.i
  let first = p.atom.first
  p.atom.first = false
  case c
  of '*':
    inc p.i
    return atomPrim(apAny)
  of '#':
    inc p.i
    let z = p.readInt
    if z < 0 or not isValidElement(z): p.fail("expected atomic number after '#'", start)
    return atomPrim(apAtomicNum, z)
  of '0'..'9':
    p.atom.first = first   # an isotope may precede [2H]
    return atomPrim(apIsotope, p.readInt)
  of '+', '-':
    inc p.i
    var n = p.readInt
    if n < 0:
      n = 1
      while p.peek == c:
        inc n
        inc p.i
    return atomPrim(apCharge, if c == '+': n else: -n)
  of '@':
    inc p.i
    p.atom.chiral = ctCCW
    if p.peek == '@':
      inc p.i
      p.atom.chiral = ctCW
    if p.peek == '?':
      p.fail("'@?' (chiral or unspecified) is not supported")
    # a chirality class; other capitals are primitives, as in [C@D3]
    let cls = p.peek & p.peek(1)
    if cls == "TH":
      p.i += 2
      case p.readInt
      of 1: p.atom.chiral = ctCCW
      of 2: p.atom.chiral = ctCW
      else: p.fail("invalid @TH chirality class", start)
    elif cls in ["AL", "SP", "SQ", "TB", "OH"]:
      p.fail("only tetrahedral chirality (@, @@) is supported", start)
    return nil
  of ':':
    inc p.i
    p.atom.atomMap = p.readInt
    if p.atom.atomMap < 0: p.fail("expected atom map number after ':'", start)
    return nil
  of '$':
    if p.peek(1) != '(': p.fail("expected '(' after '$'")
    return p.parseRecursive()
  of 'a'..'z':
    for sym in [c & p.peek(1), $c]:
      if sym in aromaticSymbols:
        p.i += sym.len
        return elementExpr(atomicNumber(sym[0].toUpperAscii & sym[1..^1]), true)
    inc p.i
    case c
    of 'a': return atomPrim(apAromatic)
    of 'h': return p.countPrim(apImplicitH)
    of 'r': return p.ringPrim(apRingSize)
    of 'v': return atomPrim(apValence, p.readIntOr(1))
    of 'x': return p.countPrim(apRingBonds)
    else: p.fail("unknown atom primitive '" & c & "'", start)
  of 'A'..'Z':
    # [H], [2H], [H+]: hydrogen atoms; otherwise H is a hydrogen count
    if c == 'H' and first and p.peek(1) notin {'0'..'9', 'a'..'z'}:
      inc p.i
      return atomPrim(apAtomicNum, 1)
    if p.peek(1).isLowerAscii:
      let z = atomicNumber(c & p.peek(1))
      if z > 0:
        p.i += 2
        return elementExpr(z, false)
    inc p.i
    case c
    of 'A': return atomPrim(apAliphatic)
    of 'D': return atomPrim(apDegree, p.readIntOr(1))
    of 'X': return atomPrim(apConnectivity, p.readIntOr(1))
    of 'H':
      let n = p.readIntOr(1)
      p.atom.hasH = n > 0
      return atomPrim(apTotalH, n)
    of 'R': return p.ringPrim(apRingCount)
    else:
      let z = atomicNumber($c)
      if z <= 0: p.fail("unknown atom primitive '" & c & "'", start)
      return elementExpr(z, false)
  else:
    p.fail("unexpected character '" & c & "' in atom expression", start)

proc parseBracketAtom(p: var Parser): (QueryAtom, bool) =
  ## The atom, and whether it has an H slot for its chirality.
  let start = p.i
  inc p.i   # '['
  p.atom = AtomState(first: true)
  let e = p.parseExpr(parseAtomPrim,
                      proc (c: char): bool = c notin {',', ';', ']', '\0'},
                      AtomPrim(kind: apAny))
  if p.peek != ']': p.fail("unterminated bracket atom", start)
  inc p.i
  (QueryAtom(expr: e, chiral: p.atom.chiral, atomMap: p.atom.atomMap),
   p.atom.chiral != ctNone and p.atom.hasH)

proc parseOrganicAtom(p: var Parser): QueryAtom =
  let c = p.peek
  let two = c & p.peek(1)
  var e: Expr[AtomPrim]
  if two in organicSubset:
    p.i += 2
    e = elementExpr(atomicNumber(two), false)
  else:
    inc p.i
    if $c in organicSubset: e = elementExpr(atomicNumber($c), false)
    elif $c in organicAromatic: e = elementExpr(atomicNumber($c.toUpperAscii), true)
    elif c == '*': e = atomPrim(apAny)
    elif c == 'a': e = atomPrim(apAromatic)
    elif c == 'A': e = atomPrim(apAliphatic)
    else: p.fail("unexpected character '" & c & "'", p.i - 1)
  QueryAtom(expr: e)

# ---------------------------------------------------------------------------
# Bond expressions

const bondChars = {'-', '=', '#', '$', ':', '~', '@', '/', '\\', '!', '&', ',', ';'}

proc parseBondPrim(p: var Parser): Expr[BondPrim] =
  let c = p.peek
  inc p.i
  case c
  of '-', '/', '\\': prim(bpSingle)   # marks: see `assignDoubleBondStereo`
  of '=': prim(bpDouble)
  of '#': prim(bpTriple)
  of '$': prim(bpQuadruple)
  of ':': prim(bpAromatic)
  of '~': prim(bpAny)
  of '@': prim(bpRing)
  else: p.fail("unexpected character '" & c & "' in bond expression", p.i - 1)

proc parseBondExpr(p: var Parser): Expr[BondPrim] =
  ## `p.i` is left after the expression.
  p.parseExpr(parseBondPrim,
              proc (c: char): bool = c in bondChars and c notin {',', ';'},
              bpAny)

# ---------------------------------------------------------------------------
# Graph

proc parseSmartsInto(p: var Parser) =
  var prev = -1
  var branchStack: seq[int]
  var pending: Expr[BondPrim] = nil
  var pendingText = ""
  var pendingDir = dNone
  var pendingPos = -1

  template clearPending() =
    pending = nil
    pendingText = ""
    pendingDir = dNone
    pendingPos = -1

  proc connect(p: var Parser; a, b: int; e: Expr[BondPrim]; dir: Dir; pos: int) =
    ## Adds bond a-b; a `/` `\` mark `dir` is as written from `a`.
    if p.q.bondBetween(a, b) >= 0:
      p.fail("atoms are already bonded", pos)
    p.q.addBond(a, b, if e.isNil: prim(bpDefault) else: e)
    p.marks.add dir
    p.bondPos.add pos

  while p.i < p.s.len:
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
      inc p.q.atoms[prev].ringDigits
      var r: RingOpen
      if p.openRings.pop(num, r):
        if r.atom == prev: p.fail("ring closure to the same atom", pos)
        # as in SMILES, both ends may give the bond, but not differently
        # (up to direction marks, which are checked below)
        if not r.bond.isNil and not pending.isNil and
           r.bondText.multiReplace(("/", "-"), ("\\", "-")) !=
           pendingText.multiReplace(("/", "-"), ("\\", "-")):
          p.fail("conflicting ring closure bond expressions", pos)
        let e = if pending.isNil: r.bond else: pending
        # each mark reads from its own atom (C/1 ... C\1 agree)
        if r.dir != dNone and pendingDir != dNone and r.dir == pendingDir:
          p.fail("conflicting ring closure bond directions", pos)
        if r.dir != dNone or pendingDir == dNone:
          p.connect(r.atom, prev, e, r.dir, pos)
        else:
          p.connect(prev, r.atom, e, pendingDir, pos)
        p.q.atoms[r.atom].nbrOrder[r.slot] = prev
        p.q.atoms[prev].nbrOrder.add r.atom
      else:
        p.openRings[num] = RingOpen(atom: prev, bond: pending, bondText: pendingText,
                                    dir: pendingDir,
                                    slot: p.q.atoms[prev].nbrOrder.len)
        p.q.atoms[prev].nbrOrder.add -2   # placeholder until the ring closes
      clearPending()
    elif c in bondChars:
      if pendingPos >= 0: p.fail("two consecutive bond expressions")
      pendingPos = p.i
      pending = p.parseBondExpr()
      pendingText = p.s[pendingPos ..< p.i]
      pendingDir = if pendingText == "/": dUp elif pendingText == "\\": dDown else: dNone
    else:
      let atomStart = p.i
      let (a, hasH) = if c == '[': p.parseBracketAtom()
                      else: (p.parseOrganicAtom(), false)
      let idx = p.q.addAtom(a)
      if prev >= 0:
        p.connect(prev, idx, pending, pendingDir,
                  if pendingPos >= 0: pendingPos else: atomStart)
        p.q.atoms[prev].nbrOrder.add idx
        p.q.atoms[idx].nbrOrder.add prev
      elif pendingPos >= 0:
        p.fail("bond symbol without a preceding atom", pendingPos)
      if hasH: p.q.atoms[idx].nbrOrder.add -1   # after the preceding atom, as in SMILES
      clearPending()
      prev = idx

  if pendingPos >= 0: p.fail("dangling bond symbol at end of input", pendingPos)
  if branchStack.len > 0: p.fail("unclosed branch '('")
  for num in p.openRings.keys:
    p.fail("unclosed ring bond " & $num)
  if p.q.numAtoms == 0: p.fail("empty SMARTS", 0)
  p.assignDoubleBondStereo()

proc assignDoubleBondStereo(p: var Parser) =
  ## Cis/trans from `/` `\` marks around bonds that are plain `=` queries.
  let isDouble = p.q.bonds.mapIt(it.expr.op == eoPrim and it.expr.prim == bpDouble)
  for st in stereoFromMarks(p.q, isDouble, p.marks):
    if st.conflict >= 0:
      p.fail("conflicting '/' '\\' marks on two substituents of the same " &
             "double-bond atom", p.bondPos[st.conflict])
    p.q.bonds[st.bond].stereo = if st.cis: bsCis else: bsTrans
    p.q.bonds[st.bond].stereoAtoms = st.refs

proc parseSmarts*(s: string): Query =
  ## Parses a SMARTS pattern into a `Query`. Raises `SmartsError` (with a
  ## character position) on syntax errors.
  var p = Parser(s: s)
  p.parseSmartsInto()
  result = p.q
  result.numRecursive = p.nRecursive
