## SELFIES (Krenn et al., 2020): a string notation for molecules in which
## every sequence of symbols decodes to a valid molecule. Compatible with
## the reference `selfies` package (version 2): `selfiesToSmiles` gives the
## same SMILES as its `decoder`, and `smilesToSelfies` the same SELFIES as
## its `encoder` for a Kekulé SMILES.
##
## A SELFIES string is a sequence of bracketed symbols:
## - atoms, with the bond to the previous atom: `[C]`, `[=O]`, `[#N]`,
##   `[/C]`, `[NH1+1]`, `[C@@H1]`, `[13CH4]`;
## - `[Branch1]`..`[Branch3]` (`[=Branch1]`, `[#Branch1]` ...): the next 1-3
##   symbols, read as base-16 digits (`indexAlphabet`), give the branch
##   length;
## - `[Ring1]`..`[Ring3]` (`[=Ring1]`, `[/\Ring1]` ...): the next 1-3
##   symbols give how far back the ring bond reaches;
## - `.` between fragments.
## Decoding tracks each atom's remaining bonding capacity (the semantic
## constraints, `bondingCapacity`) and lowers bond orders that don't fit,
## which is what makes every string valid.

import std/[strutils, tables, sets, algorithm, sequtils]
import ./[mol, elements]
import ./smiles/[parser, writer]

type
  SelfiesError* = object of ChemError
  SelfiesConstraints* = enum
    ## the reference package's preset semantic constraints
    scDefault, scOctetRule, scHypervalent

const
  indexAlphabet* = ["[C]", "[Ring1]", "[Ring2]", "[Branch1]", "[=Branch1]",
                    "[#Branch1]", "[Branch2]", "[=Branch2]", "[#Branch2]",
                    "[O]", "[N]", "[=N]", "[=C]", "[#C]", "[S]", "[P]"]
    ## the symbols standing for the digits 0-15 in branch and ring indices
  selfiesElements = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V",
    "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
    "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
    "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "Hf", "Ta", "W",
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Fl",
    "Lv", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
    "Er", "Tm", "Yb", "Lu", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"]
    ## elements the reference package accepts
  bondSymbols = [("", 1), ("=", 2), ("#", 3)]
  defaultCapacities = {"H": 1, "F": 1, "Cl": 1, "Br": 1, "I": 1,
    "B": 3, "B+1": 2, "B-1": 4, "O": 2, "O+1": 3, "O-1": 1,
    "N": 3, "N+1": 4, "N-1": 2, "C": 4, "C+1": 3, "C-1": 3,
    "P": 5, "P+1": 4, "P-1": 6, "S": 6, "S+1": 5, "S-1": 5}.toTable
  otherCapacity = 8   ## for atoms not in the table
  maxIndexSymbols = 3 ## [Ring3], [Branch3]: indices below 16^3
  maxNumberDigits = 9 ## longer isotopes or charges are rejected (int range)

func chargeKey(element: string; charge: int): string =
  if charge == 0: element
  elif charge > 0: element & "+" & $charge
  else: element & "-" & $(-charge)

func capacityOfKey(key: string; constraints: SelfiesConstraints): int =
  case constraints
  of scDefault: discard
  of scOctetRule:
    case key
    of "S": return 2
    of "S+1": return 3
    of "S-1": return 1
    of "P": return 3
    of "P+1": return 4
    of "P-1": return 2
    else: discard
  of scHypervalent:
    case key
    of "Cl", "Br", "I": return 7
    of "N": return 5
    else: discard
  defaultCapacities.getOrDefault(key, otherCapacity)

func bondingCapacity*(element: string; charge = 0;
                      constraints = scDefault): int =
  ## The most bonds (counted by order, H included) an atom may have under
  ## the semantic constraints.
  capacityOfKey(chargeKey(element, charge), constraints)

func robustAlphabet*(constraints = scDefault): seq[string] =
  ## The symbols the semantic constraints cover (the reference package's
  ## `get_semantic_robust_alphabet`), sorted: a vocabulary for generating
  ## SELFIES.
  var s = initHashSet[string]()
  for key in defaultCapacities.keys:
    let cap = capacityOfKey(key, constraints)
    for (b, order) in bondSymbols:
      if order <= cap: s.incl "[" & b & key & "]"
  for i in 1..3:
    for sym in ["Ring", "=Ring", "Branch", "=Branch", "#Branch"]:
      s.incl "[" & sym & $i & "]"
  for sym in indexAlphabet: s.incl sym
  result = toSeq(s)
  result.sort()

# ---------------------------------------------------------------------------
# The molecular graph SELFIES and SMILES are derived from: atoms in string
# order, and per atom its outgoing bonds in the order they are written (a
# ring bond appears at both of its atoms).

type
  GAtom = object
    element: string
    bare: bool          ## written without brackets (implicit H)
    isotope: int        ## -1: none
    chirality: string   ## "", "@" or "@@"
    hCount: int         ## bracket atoms only
    charge: int
  GBond = object
    dst: int            ## -1: placeholder for a ring bond not yet closed
    order: int
    stereo: char        ## '\0', '/' or '\\'
    ring: bool
  Graph = object
    atoms: seq[GAtom]
    outs: seq[seq[GBond]]
    counts: seq[int]    ## bond orders at each atom
    roots: seq[int]

func capacity(a: GAtom; c: SelfiesConstraints): int =
  bondingCapacity(a.element, a.charge, c) - (if a.bare: 0 else: a.hCount)

proc addAtom(g: var Graph; a: GAtom; root: bool): int =
  result = g.atoms.len
  if root: g.roots.add result
  g.atoms.add a
  g.outs.add @[]
  g.counts.add 0

proc insertBond(g: var Graph; src: int; b: GBond; pos: int) =
  if pos < 0 or pos == g.outs[src].len: g.outs[src].add b
  elif g.outs[src][pos].dst < 0: g.outs[src][pos] = b
  else: g.outs[src].insert(b, pos)

proc addBond(g: var Graph; src, dst, order: int; stereo: char) =
  g.insertBond(src, GBond(dst: dst, order: order, stereo: stereo), -1)
  g.counts[src] += order
  g.counts[dst] += order

proc addRingBond(g: var Graph; a, b, order: int; aStereo, bStereo: char;
                 aPos, bPos: int) =
  g.insertBond(a, GBond(dst: b, order: order, stereo: aStereo, ring: true), aPos)
  g.insertBond(b, GBond(dst: a, order: order, stereo: bStereo, ring: true), bPos)
  g.counts[a] += order
  g.counts[b] += order

func findBond(g: Graph; src, dst: int): int =
  for k, b in g.outs[src]:
    if b.dst == dst: return k
  -1

func bondChar(order: int; stereo: char): string =
  case order
  of 2: "="
  of 3: "#"
  else: (if stereo == '\0': "" else: $stereo)

func atomText(a: GAtom): string =
  ## The atom as written in SMILES or (without brackets) SELFIES.
  if a.bare: return a.element
  if a.isotope >= 0: result.add $a.isotope
  result.add a.element
  result.add a.chirality
  if a.hCount != 0: result.add "H" & $a.hCount
  elif a.isotope < 0 and a.chirality == "" and a.charge == 0 and
       a.element in organicSubset:
    result.add "H0"
  if a.charge > 0: result.add "+" & $a.charge
  elif a.charge < 0: result.add "-" & $(-a.charge)

func stereoOf(c: char): char =
  if c in {'/', '\\'}: c else: '\0'

# ---------------------------------------------------------------------------
# Decoding

iterator selfiesSymbols*(s: string): string =
  ## The symbols of a SELFIES string, with `.` between fragments (the
  ## reference package's `split_selfies`). Raises `SelfiesError` for an
  ## unclosed `[`.
  var left = s.find('[')
  while left >= 0 and left < s.len:
    let right = s.find(']', left + 1)
    if right < 0: raise newException(SelfiesError, "malformed SELFIES string, hanging '[' bracket")
    yield s[left .. right]
    left = right + 1
    if left < s.len and s[left] == '.':
      yield "."
      inc left

func bondOrder(c: char): int =
  case c
  of '=': 2
  of '#': 3
  else: 1

func scanBracketAtom(s: string; i: int; smiles: bool; atom: var GAtom): bool =
  ## The bracket atom from `s[i]` to the closing `]` (the last character):
  ## isotope, element, `@`/`@@`, H count, charge. In SELFIES the H count is
  ## `H` and one digit and the charge `+1`, `-2`...; in SMILES (`smiles`)
  ## `H` alone means one, charges may be `+`, `++` or `+2`, and an atom map
  ## (dropped) may follow.
  let stop = s.high
  var i = i
  let isoStart = i
  while i < stop and s[i] in Digits: inc i
  let iso = s[isoStart ..< i]
  if iso.len > maxNumberDigits or i >= stop or s[i] notin {'A'..'Z'}: return false
  var element = $s[i]
  inc i
  if i < stop and s[i] in {'a'..'z'}:
    element.add s[i]
    inc i
  var chir = ""
  while i < stop and s[i] == '@' and chir.len < 2:
    chir.add '@'
    inc i
  var h = 0
  if i < stop and s[i] == 'H':
    inc i
    if i < stop and s[i] in Digits:
      h = ord(s[i]) - ord('0')
      inc i
    elif smiles: h = 1
    else: return false
  var charge = 0
  if i < stop and s[i] in {'+', '-'}:
    let c = s[i]
    let sign = if c == '+': 1 else: -1
    inc i
    let start = i
    if smiles:
      var n = 1
      while i < stop and s[i] == c:
        inc n
        inc i
      if n == 1 and i < stop and s[i] in Digits:
        while i < stop and s[i] in Digits: inc i
        if i - start > maxNumberDigits: return false
        n = parseInt(s[start ..< i])
      charge = sign * n
    else:
      while i < stop and s[i] in {'1'..'9'}: inc i
      if i == start or i - start > maxNumberDigits: return false
      charge = sign * parseInt(s[start ..< i])
  if smiles and i < stop and s[i] == ':':   # atom map: not kept
    inc i
    if i >= stop or s[i] notin Digits: return false
    while i < stop and s[i] in Digits: inc i
  if i != stop or element notin selfiesElements: return false
  atom = GAtom(element: element, isotope: (if iso.len > 0: parseInt(iso) else: -1),
               chirality: chir, hCount: h, charge: charge)
  true

func parseAtomSymbol(sym: string; atom: var GAtom; order: var int;
                     stereo: var char): bool =
  ## `[` bond? then an organic-subset element, or a bracket atom
  if sym.len < 3 or sym[0] != '[' or sym[^1] != ']': return false
  let bondC = if sym[1] in {'=', '#', '/', '\\'}: sym[1] else: '\0'
  let i = if bondC == '\0': 1 else: 2
  order = bondOrder(bondC)
  stereo = stereoOf(bondC)
  let inner = sym[i ..< sym.high]
  if inner in organicSubset:
    atom = GAtom(element: inner, bare: true, isotope: -1)
    return true
  scanBracketAtom(sym, i, smiles = false, atom)

func buildBranchSymbols(): Table[string, (int, int)] =
  for (b, order) in bondSymbols:
    for l in 1..3: result["[" & b & "Branch" & $l & "]"] = (order, l)

func buildRingSymbols(): Table[string, (int, int, (char, char))] =
  for l in 1..3:
    for (b, order) in bondSymbols:
      result["[" & b & "Ring" & $l & "]"] = (order, l, ('\0', '\0'))
    for lc in ['-', '/', '\\']:
      for rc in ['-', '/', '\\']:
        if lc != '-' or rc != '-':
          result["[" & lc & rc & "Ring" & $l & "]"] = (1, l, (stereoOf(lc), stereoOf(rc)))

const
  branchSymbols = buildBranchSymbols()   ## symbol -> (bond order, index length)
  ringSymbols = buildRingSymbols()       ## ... and the stereo marks at each end

type
  PendingRing = tuple[a, b, order: int; aStereo, bStereo: char]
  Decoder = object
    g: Graph
    toks: seq[string]
    pos: int
    rings: seq[PendingRing]
    c: SelfiesConstraints
    text: string

proc fail(d: Decoder; sym: string) {.noreturn.} =
  raise newException(SelfiesError, "invalid symbol '" & sym & "' in SELFIES " & d.text)

proc readIndex(d: var Decoder; n: int): int =
  ## The next `n` symbols as base-16 digits (missing or unknown ones: 0).
  for _ in 1..n:
    var digit = 0
    if d.pos < d.toks.len:
      digit = max(0, indexAlphabet.find(d.toks[d.pos]))
      inc d.pos
    result = result * indexAlphabet.len + digit

type Frame = object
  ## a derivation in progress: the main chain of a fragment, or a branch
  budget: int     ## how many symbols it may read
  used: int       ## how many it has read
  state: int      ## bonds still available to `prev`
  prev: int       ## the atom the next symbol attaches to
  alive: bool

proc deriveFragment(d: var Decoder) =
  ## Reads one fragment's symbols. A branch reads its own number of
  ## symbols, as a frame on an explicit stack (branches can nest
  ## arbitrarily deep). Each frame tracks its atom's remaining bonds:
  ## bond orders are lowered to fit, and a full atom ends the frame.
  var stack = @[Frame(budget: high(int), prev: -1, alive: true)]
  while stack.len > 0:
    if not (stack[^1].alive and stack[^1].used < stack[^1].budget and d.pos < d.toks.len):
      # done: the rest of its budget is skipped, and counts for the parent
      let f = stack.pop()
      var used = f.used
      while used < f.budget and d.pos < d.toks.len:
        inc d.pos
        inc used
      if stack.len > 0: stack[^1].used += used
      continue
    let sym = d.toks[d.pos]
    inc d.pos
    inc stack[^1].used
    let state = stack[^1].state
    let prev = stack[^1].prev
    # as the reference tells symbols apart: by their 4th- and 3rd-last
    # characters ("[=Branch1]", "[Ring2]")
    let tag = if sym.len >= 4: (sym[^4], sym[^3]) else: ('\0', '\0')
    if tag == ('c', 'h'):
      let (btype, n) = branchSymbols.getOrDefault(sym)
      if n == 0: d.fail(sym)
      if state > 1:
        let branchState = min(state - 1, btype)
        stack[^1].state = state - branchState
        let q = d.readIndex(n)
        stack[^1].used += n
        stack.add Frame(budget: q + 1, state: branchState, prev: prev, alive: true)
    elif tag == ('n', 'g'):
      let (rtype, n, stereo) = ringSymbols.getOrDefault(sym)
      if n == 0: d.fail(sym)
      if state != 0:
        let order = min(rtype, state)
        stack[^1].state = state - order
        stack[^1].alive = state - order != 0
        let q = d.readIndex(n)
        stack[^1].used += n
        d.rings.add (max(0, prev - (q + 1)), prev, order, stereo[0], stereo[1])
    elif "eps" in sym:
      if state != 0: stack[^1].alive = false
    else:
      var atom: GAtom
      var order: int
      var stereo: char
      if not parseAtomSymbol(sym, atom, order, stereo): d.fail(sym)
      let cap = atom.capacity(d.c)
      if cap < 0: d.fail(sym)   # more H than the atom can have
      order = min([order, state, cap])   # 0 for a fragment's first atom
      stack[^1].state = cap - order
      stack[^1].alive = cap - order != 0
      if order == 0:
        if state == 0: stack[^1].prev = d.g.addAtom(atom, root = true)
      else:
        let i = d.g.addAtom(atom, root = false)
        d.g.addBond(prev, i, order, stereo)
        stack[^1].prev = i

proc formRings(d: var Decoder) =
  ## Ring bonds, once all atoms are known: limited by both atoms' free
  ## capacity; a second ring bond between two atoms raises the bond order.
  var made = newSeq[int](d.g.atoms.len)
  for (a, b, order0, aStereo, bStereo) in d.rings:
    if a == b: continue
    let aFree = d.g.atoms[a].capacity(d.c) - d.g.counts[a]
    let bFree = d.g.atoms[b].capacity(d.c) - d.g.counts[b]
    if aFree <= 0 or bFree <= 0: continue
    let order = min([order0, aFree, bFree])
    let k = d.g.findBond(a, b)
    if k >= 0:
      let old = d.g.outs[a][k].order
      let newOrder = min(order + old, 3)
      d.g.outs[a][k].order = newOrder
      if d.g.outs[a][k].ring: d.g.outs[b][d.g.findBond(b, a)].order = newOrder
      d.g.counts[a] += newOrder - old
      d.g.counts[b] += newOrder - old
    else:
      d.g.addRingBond(a, b, order, aStereo, bStereo, made[a], made[b])
      inc made[a]
      inc made[b]

proc toSmilesText(g: Graph): string =
  ## The SMILES for the graph in its order, as the reference decoder writes
  ## it. Ring-bond numbers are not reused (as there) until they would pass
  ## 99, where the reference writes invalid `%100`; from then on the lowest
  ## free number is used.
  var ringNum = initTable[(int, int), int]()
  var openNums = initHashSet[int]()
  var parts: seq[string]
  for root in g.roots:
    var s = ""
    var stack = @[(root, 0, false)]   # atom, next out-bond, close a branch
    while stack.len > 0:
      let (cur, k, closing) = stack[^1]
      let outs = g.outs[cur]
      if k == 0: s.add (if g.atoms[cur].bare: atomText(g.atoms[cur])
                        else: "[" & atomText(g.atoms[cur]) & "]")
      if k < outs.len:
        let b = outs[k]
        stack[^1] = (cur, k + 1, closing)
        if b.ring:
          s.add bondChar(b.order, b.stereo)
          let ends = (min(cur, b.dst), max(cur, b.dst))
          var num = ringNum.getOrDefault(ends, 0)
          if num == 0:
            num = ringNum.len + 1
            if num > 99:
              num = 1
              while num in openNums: inc num
            ringNum[ends] = num
            openNums.incl num
          else:
            openNums.excl num
          if num >= 10: s.add '%'
          s.add $num
        else:
          let branch = k < outs.len - 1
          if branch: s.add '('
          s.add bondChar(b.order, b.stereo)
          stack.add (b.dst, 0, branch)
      else:
        discard stack.pop()
        if closing: s.add ')'
    parts.add s
  parts.join(".")

proc selfiesToSmiles*(s: string; constraints = scDefault): string =
  ## Decodes SELFIES to SMILES, as the reference package's `decoder`:
  ## Kekulé form, atoms in SELFIES order. Raises `SelfiesError` for a
  ## symbol that isn't valid SELFIES (any sequence of valid symbols
  ## decodes).
  var d = Decoder(c: constraints, text: s)
  for fragment in s.split('.'):
    d.toks.setLen 0
    d.pos = 0
    for sym in selfiesSymbols(fragment):
      if sym != "[nop]": d.toks.add sym
    d.deriveFragment()
  d.formRings()
  toSmilesText(d.g)

proc parseSelfies*(s: string; sanitize = true; removeHs = true;
                   constraints = scDefault): Mol =
  ## Reads a molecule from SELFIES (see `selfiesToSmiles` and
  ## `parseSmiles`).
  parseSmiles(selfiesToSmiles(s, constraints), sanitize, removeHs)

# ---------------------------------------------------------------------------
# Encoding

proc fail(msg, smiles: string) {.noreturn.} =
  raise newException(SelfiesError, msg & " in SMILES " & smiles)

proc parseSmilesAtom(tok: string; a: var GAtom): bool =
  ## An atom token of a Kekulé SMILES: organic subset or bracket atom.
  if tok in organicSubset:
    a = GAtom(element: tok, bare: true, isotope: -1)
    return true
  tok.len >= 3 and tok[0] == '[' and tok[^1] == ']' and
    scanBracketAtom(tok, 1, smiles = true, a)

type AromaticSmiles = object of SelfiesError
  ## the SMILES has aromatic atoms or bonds: kekulize it first

func isAromaticAtom(tok: string): bool =
  ## `c`, `[nH]`, `[13c]`...
  var i = if tok[0] == '[': 1 else: 0
  while i < tok.len and tok[i] in Digits: inc i
  i < tok.len and tok[i] in {'a'..'z'}

proc kekuleSmilesGraph(smiles: string): Graph =
  ## The graph of a Kekulé SMILES in written order (the reference
  ## package's `smiles_to_mol`).
  var i = 0
  while i < smiles.len:
    # one fragment
    var prevStack = @[-1]
    var branches = 0
    var rings = initTable[string, (char, int, int)]()   # bond char, atom, slot
    var chainStart = true
    while i < smiles.len:
      if smiles[i] == '.':
        inc i
        break
      var bondC = '\0'
      if smiles[i] in {'-', '/', '\\', '=', '#'}:
        bondC = smiles[i]
        inc i
        if i == smiles.len: fail("hanging bond", smiles)
      elif smiles[i] == ':':
        raise newException(AromaticSmiles, "aromatic bond")
      elif smiles[i] == '$':
        fail("quadruple bonds are not supported", smiles)
      let c = smiles[i]
      let prev = prevStack[^1]
      if c in Letters or c == '[':
        var tok: string
        if c == '[':
          let r = smiles.find(']', i + 1)
          if r < 0: fail("hanging '['", smiles)
          tok = smiles[i .. r]
        elif smiles.continuesWith("Br", i) or smiles.continuesWith("Cl", i):
          tok = smiles[i .. i + 1]
        else:
          tok = $c
        i += tok.len
        var a: GAtom
        if not parseSmilesAtom(tok, a):
          if isAromaticAtom(tok): raise newException(AromaticSmiles, "aromatic atom " & tok)
          fail("invalid or unsupported atom " & tok & " (SELFIES has no wildcard atoms)", smiles)
        let idx = result.addAtom(a, root = prev < 0)
        if prev >= 0: result.addBond(prev, idx, bondOrder(bondC), stereoOf(bondC))
        prevStack[^1] = idx
        chainStart = false
      elif chainStart:
        fail("chain begins with '" & c & "'", smiles)
      elif c in {'(', ')'}:
        if bondC != '\0': fail("hanging bond", smiles)
        inc i
        if c == '(':
          inc branches
          prevStack.add prev
          chainStart = true
        else:
          if branches == 0: fail("hanging ')'", smiles)
          dec branches
          discard prevStack.pop()
      elif c in Digits or c == '%':
        let num = if c == '%': smiles[i ..< min(i + 3, smiles.len)] else: $c
        if c == '%' and (num.len != 3 or num[1] notin Digits or num[2] notin Digits):
          fail("invalid ring number", smiles)
        i += num.len
        if num notin rings:
          result.outs[prev].add GBond(dst: -1)   # its slot
          rings[num] = (bondC, prev, result.outs[prev].high)
        else:
          let (lc, latom, lpos) = rings[num]
          rings.del num
          if latom == prev: fail("ring bond from an atom to itself", smiles)
          if result.findBond(latom, prev) >= 0 or result.findBond(prev, latom) >= 0:
            fail("ring bond between bonded atoms", smiles)
          var (x, y) = (lc, bondC)
          if x == '\0': swap(x, y)
          if not (x == y or y == '\0' or (x in {'/', '\\'} and y in {'/', '\\'})):
            fail("mismatched ring bonds", smiles)
          result.addRingBond(latom, prev, max(bondOrder(lc), bondOrder(bondC)),
                             stereoOf(lc), stereoOf(bondC), lpos, -1)
      else:
        fail("unexpected '" & c & "'", smiles)
    if branches > 0: fail("hanging '('", smiles)
    if rings.len > 0: fail("unclosed ring bond", smiles)

func indexLength(index: int): int =
  ## How many symbols `index` takes in base 16.
  result = 1
  var x = index div indexAlphabet.len
  while x > 0:
    inc result
    x = x div indexAlphabet.len

proc addIndex(s: var string; index: int) =
  ## `index` in base 16, most significant digit first, as SELFIES symbols.
  var place = 1
  for _ in 2..indexLength(index): place *= indexAlphabet.len
  while place > 0:
    s.add indexAlphabet[index div place mod indexAlphabet.len]
    place = place div indexAlphabet.len

proc shouldInvertChirality(g: Graph; atom: int): bool =
  ## The decoder lists an atom's ring bonds first, in the order they are
  ## closed; whether that reorders its neighbours by an odd permutation.
  var closing, opening, other: seq[int]
  for k, b in g.outs[atom]:
    if not b.ring: other.add k
    elif atom < b.dst: opening.add k
    else: closing.add k
  opening.sort(proc (x, y: int): int = cmp(g.outs[atom][x].dst, g.outs[atom][y].dst))
  let perm = closing & opening & other
  permutationParity(toSeq(0..<perm.len), perm) == 1

proc fragmentLengths(g: Graph; smiles: string): seq[int] =
  ## Per atom, how many symbols the fragment written from it has (the atom,
  ## its rings and branches, and the chain after it). A child atom always
  ## comes after its parent, so one pass from the end finds them all.
  result = newSeq[int](g.atoms.len)
  for a in countdown(g.atoms.high, 0):
    var n = 1
    for k, b in g.outs[a]:
      if b.ring:
        if a > b.dst:   # written at the ring's later atom
          let digits = indexLength(a - b.dst - 1)
          if digits > maxIndexSymbols: fail("ring bond too long for SELFIES", smiles)
          n += 1 + digits
      elif k == g.outs[a].high:
        n += result[b.dst]
      else:
        let digits = indexLength(result[b.dst] - 1)
        if digits > maxIndexSymbols: fail("branch too long for SELFIES", smiles)
        n += 1 + digits + result[b.dst]
    result[a] = n

proc writeFragment(g: Graph; root: int; lengths: seq[int]; s: var string) =
  ## Writes the fragment from `root`: each atom, then its out-bonds in
  ## order (ring bonds at their later atom; branches, with their length
  ## first; the last bond continues the chain).
  var stack = @[(root, 0, GBond(order: 1))]   # atom, next out-bond, bond in
  while stack.len > 0:
    let (cur, k, incoming) = stack[^1]
    if k == 0:
      s.add '['
      s.add bondChar(incoming.order, incoming.stereo)
      s.add atomText(g.atoms[cur])
      s.add ']'
    if k == g.outs[cur].len:
      discard stack.pop()
      continue
    let b = g.outs[cur][k]
    stack[^1][1] = k + 1
    if b.ring:
      if cur > b.dst:
        let back = g.outs[b.dst][g.findBond(b.dst, cur)]
        s.add '['
        if back.order != 1 or (back.stereo == '\0' and b.stereo == '\0'):
          s.add bondChar(back.order, '\0')
        else:
          s.add (if back.stereo == '\0': '-' else: back.stereo)
          s.add (if b.stereo == '\0': '-' else: b.stereo)
        s.add "Ring" & $indexLength(cur - b.dst - 1) & "]"
        s.addIndex(cur - b.dst - 1)
    elif k == g.outs[cur].high:
      stack[^1] = (b.dst, 0, b)   # the chain goes on
    else:
      let n = lengths[b.dst]
      s.add "[" & bondChar(b.order, '\0') & "Branch" & $indexLength(n - 1) & "]"
      s.addIndex(n - 1)
      stack.add (b.dst, 0, b)

proc kekuleSmilesToSelfies(smiles: string; strict: bool;
                           constraints: SelfiesConstraints): string =
  var g = kekuleSmilesGraph(smiles)
  if strict:
    for i, a in g.atoms:
      if g.counts[i] > a.capacity(constraints):
        fail("atom " & atomText(a) & " has " & $g.counts[i] & " bonds, more than the " &
             $a.capacity(constraints) & " SELFIES allows", smiles)
  for i in 0..<g.atoms.len:
    if g.atoms[i].chirality != "" and shouldInvertChirality(g, i):
      g.atoms[i].chirality = if g.atoms[i].chirality == "@": "@@" else: "@"
  let lengths = fragmentLengths(g, smiles)
  for k, root in g.roots:
    if k > 0: result.add '.'
    writeFragment(g, root, lengths, result)

proc smilesToSelfies*(smiles: string; strict = true;
                      constraints = scDefault): string =
  ## Encodes SMILES as SELFIES, keeping the SMILES atom order, as the
  ## reference package's `encoder` (whose result it equals for a Kekulé
  ## SMILES). Aromatic SMILES is first rewritten in Kekulé form, keeping
  ## its atoms (explicit H included) in order, though bracket atoms the
  ## writer doesn't need may lose their brackets (`[CH3]` as `C`). With
  ## `strict`, raises `SelfiesError` for an atom with more
  ## bonds than the semantic constraints allow (it would not decode back
  ## to the same molecule); also for wildcard atoms, quadruple bonds and
  ## ring bonds or branches too long for SELFIES.
  if smiles.len == 0: raise newException(SelfiesError, "empty SMILES")
  try:
    kekuleSmilesToSelfies(smiles, strict, constraints)
  except AromaticSmiles:
    let kekule = parseSmiles(smiles, removeHs = false).toSmiles(canonical = false, kekule = true)
    kekuleSmilesToSelfies(kekule, strict, constraints)

proc toSelfies*(m: Mol; canonical = true; strict = true;
                constraints = scDefault): string =
  ## `m` as SELFIES: from its canonical Kekulé SMILES (so the same molecule
  ## always gives the same SELFIES), or in atom order. See `smilesToSelfies`.
  kekuleSmilesToSelfies(m.toSmiles(canonical, kekule = true), strict, constraints)
