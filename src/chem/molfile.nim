## MDL molfiles (V2000 and V3000): reading and writing.
##
## Reading takes atoms, bonds, charges, isotopes, atom maps, radicals and
## explicit valences, and stereo from the coordinates (3D positions or 2D
## wedge/hash bonds) or, for files without coordinates, from atom parities.
## Writing gives a Kekulé structure with the molecule's coordinates, if it
## has any (choosing wedges for its chiral atoms in 2D), and atom parities.
## Without coordinates, double-bond cis/trans can't be written, and
## tetrahedral stereo is only in the parities, which RDKit ignores.
## Query features (atom lists, query bond types, S-groups) are not read.

import std/[strutils, sequtils, tables, math, parseutils]
import ./[mol, elements, valence, aromaticity, geometry, sanitization, stereo, canon]

type
  MolfileError* = object of ChemError
    line*: int   ## 1-based line number in the block (0 if unknown)

  RawAtom = object
    atom: Atom
    pos: Point3
    parity: int
    valence: int    ## total valence given in the file, or -1
    radical: int    ## MDL code: 1 singlet, 2 doublet, 3 triplet
    symbolIsotope: int  ## from a D or T symbol (kept when `M  ISO` resets)

  RawBond = object
    a, b: int
    order: BondOrder
    wedge: BondWedge
    unknown: bool   ## "either" stereo (V2000 3/4, V3000 CFG=2)

  Ctab = object
    atoms: seq[RawAtom]
    bonds: seq[RawBond]

proc fail(msg: string; line: int) {.noreturn.} =
  var e = newException(MolfileError, msg & (if line > 0: " (line " & $line & ")" else: ""))
  e.line = line
  raise e

func field(s: string; a, b: int): string =
  ## Columns a..b (0-based, inclusive) of a fixed-width line, stripped.
  if a >= s.len: "" else: s[a .. min(b, s.high)].strip

func trimmed(s: string; a, b: int): (int, int) =
  ## Columns a..b of a fixed-width line without surrounding spaces, as
  ## (first, last) indices; empty when last < first.
  var p = a
  var q = min(b, s.high)
  while p <= q and s[p] == ' ': inc p
  while q >= p and s[q] == ' ': dec q
  (p, q)

proc intField(s: string; a, b, line: int): int =
  ## Columns a..b as an integer (0 if blank), parsed in place.
  let (p, q) = trimmed(s, a, b)
  if p <= q and parseInt(s.toOpenArray(p, q), result) != q - p + 1:
    fail("expected a number in columns " & $(a + 1) & "-" & $(b + 1), line)

proc floatField(s: string; a, b, line: int): float =
  let (p, q) = trimmed(s, a, b)
  if p <= q and parseFloat(s.toOpenArray(p, q), result) != q - p + 1:
    fail("expected a coordinate in columns " & $(a + 1) & "-" & $(b + 1), line)

proc elementFor(sym: string; line: int): (int, int) =
  ## (atomic number, isotope) for an atom symbol; pseudo-atoms are 0.
  case sym
  of "D": (1, 2)
  of "T": (1, 3)
  of "*", "R", "R#": (0, 0)   # dummy atom, R group
  of "A", "Q", "X", "M", "L", "AH", "QH", "XH", "MH", "LP":
    fail("query atom '" & sym & "' is not supported", line)
  else:
    let z = atomicNumber(sym)
    if z < 0: fail("unknown atom symbol '" & sym & "'", line)
    (z, 0)

proc bondOrderFor(code, line: int): BondOrder =
  case code
  of 1: boSingle
  of 2: boDouble
  of 3: boTriple
  of 4: boAromatic
  else: fail("unsupported bond type " & $code & " (query or coordination bond)", line)

# ---------------------------------------------------------------------------
# V2000

const v2000Charges = [0, 3, 2, 1, 0, -1, -2, -3]   # by atom-block charge code

proc readV2000(lines: openArray[string]; first: int; ct: var Ctab): int =
  ## Reads the counts line at `first` and what follows up to `M  END`;
  ## returns the index after `M  END`.
  let counts = lines[first]
  let na = intField(counts, 0, 2, first + 1)
  let nb = intField(counts, 3, 5, first + 1)
  var i = first + 1
  for k in 0..<na:
    if i >= lines.len: fail("missing atom lines", i + 1)
    let s = lines[i]
    var ra: RawAtom
    ra.pos = (floatField(s, 0, 9, i + 1), floatField(s, 10, 19, i + 1),
              floatField(s, 20, 29, i + 1))
    let (z, iso) = elementFor(field(s, 31, 33), i + 1)
    ra.atom = initAtom(z, isotope = iso)
    ra.symbolIsotope = iso
    let massDiff = intField(s, 34, 35, i + 1)
    if massDiff != 0 and z > 0:
      ra.atom.isotope = int(round(avgMass(z))) + massDiff
    let chargeCode = intField(s, 36, 38, i + 1)
    if chargeCode notin 0..7: fail("invalid charge code " & $chargeCode, i + 1)
    ra.atom.charge = v2000Charges[chargeCode]
    if chargeCode == 4: ra.radical = 2
    ra.parity = intField(s, 39, 41, i + 1)
    ra.valence = case intField(s, 48, 50, i + 1)
      of 0: -1
      of 15: 0
      else: intField(s, 48, 50, i + 1)
    ra.atom.atomMap = intField(s, 60, 62, i + 1)
    ct.atoms.add ra
    inc i
  for k in 0..<nb:
    if i >= lines.len: fail("missing bond lines", i + 1)
    let s = lines[i]
    let a = intField(s, 0, 2, i + 1) - 1
    let b = intField(s, 3, 5, i + 1) - 1
    if a notin 0..<na or b notin 0..<na:
      fail("bond to a non-existent atom", i + 1)
    var rb = RawBond(a: a, b: b, order: bondOrderFor(intField(s, 6, 8, i + 1), i + 1))
    case intField(s, 9, 11, i + 1)
    of 1: rb.wedge = bwWedge
    of 6: rb.wedge = bwHash
    of 3, 4: rb.unknown = true
    else: discard
    ct.bonds.add rb
    inc i
  # properties block; CHG/RAD and ISO lines replace the atom-block values
  var chargesReset, isotopesReset = false
  while i < lines.len:
    let s = lines[i]
    inc i
    if s.startsWith("M  END"): return i
    if s.startsWith("A  ") or s.startsWith("G  "):
      inc i   # alias / group: its text is on the next line
    elif s.startsWith("M  ALS"):
      fail("atom lists are not supported", i)
    elif s.startsWith("S  SKP"):
      i += intField(s, 6, 8, i)
    elif s.startsWith("M  CHG") or s.startsWith("M  RAD") or s.startsWith("M  ISO"):
      let key = s[3..5]
      if key in ["CHG", "RAD"] and not chargesReset:
        for ra in ct.atoms.mitems:
          ra.atom.charge = 0
          ra.radical = 0
        chargesReset = true
      if key == "ISO" and not isotopesReset:
        for ra in ct.atoms.mitems: ra.atom.isotope = ra.symbolIsotope
        isotopesReset = true
      let n = intField(s, 6, 8, i)
      for k in 0..<n:
        let atomIdx = intField(s, 9 + 8 * k, 12 + 8 * k, i) - 1
        let value = intField(s, 13 + 8 * k, 16 + 8 * k, i)
        if atomIdx notin 0..<na: fail("property for a non-existent atom", i)
        case key
        of "CHG": ct.atoms[atomIdx].atom.charge = value
        of "RAD": ct.atoms[atomIdx].radical = value
        else: ct.atoms[atomIdx].atom.isotope = value
  fail("missing 'M  END'", i)

# ---------------------------------------------------------------------------
# V3000

proc v30Lines(lines: openArray[string]; i: var int): seq[(string, int)] =
  ## Logical `M  V30` lines (continuations joined) up to `M  END`, with
  ## their line numbers; leaves `i` after `M  END`.
  var cur = ""
  var curLine = 0
  while i < lines.len:
    let s = lines[i]
    inc i
    if s.startsWith("M  END"):
      if cur.len > 0: fail("unterminated V3000 continuation line", curLine)
      return
    if not s.startsWith("M  V30 "): continue
    if cur.len == 0: curLine = i
    var body = s[7..^1]
    if body.endsWith("-"):
      cur &= body[0 ..^ 2]
    else:
      result.add (cur & body, curLine)
      cur = ""
  fail("missing 'M  END'", i)

func v30Tokens(s: string; toks: var seq[string]) =
  ## Splits on spaces, keeping quoted strings and parenthesised lists whole.
  toks.setLen 0
  var start = -1
  var depth = 0
  var quoted = false
  for k, c in s:
    if c == ' ' and depth == 0 and not quoted:
      if start >= 0: toks.add s.substr(start, k - 1)
      start = -1
      continue
    if start < 0: start = k
    if c == '"': quoted = not quoted
    elif c == '(' and not quoted: inc depth
    elif c == ')' and not quoted: dec depth
  if start >= 0: toks.add s.substr(start)

proc intToken(t: string; line: int): int =
  try: parseInt(t)
  except ValueError: fail("expected a number, got '" & t & "'", line)

iterator v30Properties(toks: seq[string]; first, line: int): (string, int) =
  ## Numeric `KEY=value` tokens from `first` on (others, such as lists, are
  ## skipped).
  for k in first..toks.high:
    let t = toks[k]
    let eq = t.find('=')
    if eq > 0 and t.len > eq + 1 and t[eq + 1] != '(':
      yield (t[0 ..< eq], intToken(t[eq + 1 .. ^1], line))

proc readV3000(lines: openArray[string]; first: int; ct: var Ctab): int =
  var i = first + 1
  let logical = v30Lines(lines, i)
  result = i
  var section = ""
  var index: Table[int, int]   # file atom index -> position
  var skipDepth = 0            # inside a section we don't read
  var toks: seq[string]
  for (text, line) in logical:
    v30Tokens(text, toks)
    if toks.len == 0: continue
    if toks[0] == "BEGIN":
      if toks.len > 1 and toks[1] in ["ATOM", "BOND", "CTAB"] and skipDepth == 0:
        section = toks[1]
      else:
        inc skipDepth
      continue
    if toks[0] == "END":
      if skipDepth > 0: dec skipDepth
      else: section = ""
      continue
    if skipDepth > 0 or toks[0] == "COUNTS": continue
    template parseIntTok(t: string): int = intToken(t, line)
    case section
    of "ATOM":
      if toks.len < 6: fail("short V3000 atom line", line)
      var sym = toks[1]
      if sym.startsWith("\"") or sym.startsWith("[") or sym.startsWith("NOT"):
        fail("atom lists and query atoms are not supported", line)
      let (z, iso) = elementFor(sym, line)
      var ra = RawAtom(atom: initAtom(z, isotope = iso), valence: -1)
      try:
        ra.pos = (parseFloat(toks[2]), parseFloat(toks[3]), parseFloat(toks[4]))
      except ValueError: fail("expected coordinates", line)
      ra.atom.atomMap = parseIntTok(toks[5])
      for (key, v) in v30Properties(toks, 6, line):
        case key
        of "CHG": ra.atom.charge = v
        of "MASS": ra.atom.isotope = v
        of "CFG": ra.parity = v
        of "RAD": ra.radical = v
        of "VAL": ra.valence = if v == 0: -1 elif v == -1: 0 else: v   # 0: default; -1: zero
        else: discard
      index[parseIntTok(toks[0])] = ct.atoms.len
      ct.atoms.add ra
    of "BOND":
      if toks.len < 4: fail("short V3000 bond line", line)
      let a = index.getOrDefault(parseIntTok(toks[2]), -1)
      let b = index.getOrDefault(parseIntTok(toks[3]), -1)
      if a < 0 or b < 0: fail("bond to a non-existent atom", line)
      var rb = RawBond(a: a, b: b, order: bondOrderFor(parseIntTok(toks[1]), line))
      for (key, v) in v30Properties(toks, 4, line):
        if key == "CFG":
          case v
          of 1: rb.wedge = bwWedge
          of 3: rb.wedge = bwHash
          of 2: rb.unknown = true
          else: discard
      ct.bonds.add rb
    else: discard

# ---------------------------------------------------------------------------
# Building the molecule

proc toMol(ct: Ctab; sanitize: bool): Mol =
  var m: Mol
  for ra in ct.atoms: m.addAtom ra.atom
  for rb in ct.bonds:
    try: m.addBond(rb.a, rb.b, rb.order)
    except ChemError as e: fail(e.msg, 0)
    if rb.order == boAromatic:
      m.atoms[rb.a].aromatic = true
      m.atoms[rb.b].aromatic = true
  # hydrogens: an explicit valence, or fewer H for radicals
  for i, ra in ct.atoms:
    if ra.valence >= 0:
      m.atoms[i].noImplicit = true
      m.atoms[i].explicitH = max(0, ra.valence - explicitValence(m, i))
    elif ra.radical != 0:
      let unpaired = if ra.radical == 2: 1 else: 2
      m.atoms[i].explicitH = max(0, implicitHCount(m, i) - unpaired)
      m.atoms[i].noImplicit = true
  # stereo: from coordinates if there are any, else from parities
  if ct.atoms.anyIt(it.pos != (0.0, 0.0, 0.0)):
    m.coords = ct.atoms.mapIt(it.pos)
    stereoFromCoords(m, ct.bonds.mapIt(it.wedge), ct.bonds.mapIt(it.unknown))
  else:
    for i, ra in ct.atoms:
      m.atoms[i].chiral = chiralFromParity(m, i, ra.parity)
  # geometry gives every possible centre a tag: keep only real stereo,
  # even unsanitized
  if sanitize: sanitization.sanitize(m)
  else: assignStereo(m)
  m

proc readMolLines*(lines: openArray[string]; first = 0;
                   sanitize = true): (Mol, int) =
  ## Reads the molfile starting at `lines[first]`; returns the molecule and
  ## the index of the line after `M  END`. See `parseMolBlock`.
  if lines.len - first < 4: fail("molfile too short", lines.len)
  var ct: Ctab
  let version = field(lines[first + 3], 33, 38)
  let next = case version
    of "V3000": readV3000(lines, first + 3, ct)
    of "V2000", "": readV2000(lines, first + 3, ct)
    else: fail("unknown molfile version '" & version & "'", first + 4)
  (toMol(ct, sanitize), next)

proc parseMolBlock*(s: string; sanitize = true): Mol =
  ## Parses a molfile (V2000 or V3000). Raises `MolfileError`, and with
  ## `sanitize` (default) `KekulizeError` or `ValenceError` like
  ## `parseSmiles`. The coordinates are kept in `coords` (empty if they
  ## are all zero).
  readMolLines(s.splitLines, 0, sanitize)[0]

# ---------------------------------------------------------------------------
# Writing

type
  BondCfg = enum
    bcNone, bcWedge, bcHash, bcEither   ## either: unknown double-bond E/Z

  Prepared = object
    ## a Kekulé copy of the molecule, and per bond how to write it
    m: Mol
    ends: seq[(int, int)]   ## atoms, the narrow end of a wedge first
    cfg: seq[BondCfg]

func writtenValence(m: Mol; i: int): int =
  ## The total valence to write when a reader's default implicit H would be
  ## wrong, else -1.
  if defaultImplicitHCount(m, i) == totalHCount(m, i): -1 else: totalValence(m, i)

proc unspecifiedStereoBonds(m: Mol): seq[bool] =
  ## Double bonds without cis/trans that could have it (see
  ## `stereo.assignStereo`): written as "either" so that readers don't take
  ## E/Z from the drawing.
  var trial = m
  result = newSeq[bool](m.numBonds)
  for bi, bd in m.bonds:
    if bd.order != boDouble or bd.stereo != bsNone: continue
    let ra = substituents(m, bd.a, bd.b)
    let rb = substituents(m, bd.b, bd.a)
    if ra.len > 0 and rb.len > 0:
      trial.setDoubleBondStereo(bi, ra[0], rb[0], cis = true)
  assignStereo(trial)
  for bi in 0..<m.numBonds:
    result[bi] = m.bonds[bi].stereo == bsNone and trial.bonds[bi].stereo != bsNone

proc prepare(m: Mol; kekule: bool): Prepared =
  var k = m
  if kekule:
    try: kekulize(k)
    except KekulizeError: k = m   # write aromatic bonds instead
  result.m = k
  let wedges = chooseWedges(k)
  let either = if k.coords.len > 0: unspecifiedStereoBonds(k)
               else: newSeq[bool](k.numBonds)
  for bi, bd in k.bonds:
    let (narrow, w) = wedges[bi]
    result.ends.add(if narrow == bd.b: (bd.b, bd.a) else: (bd.a, bd.b))
    result.cfg.add(case w
      of bwWedge: bcWedge
      of bwHash: bcHash
      of bwNone: (if either[bi]: bcEither else: bcNone))

func dimension(m: Mol): string =
  if m.coords.len == 0: "" elif m.is3D: "3D" else: "2D"

func bondTypeCode(o: BondOrder): int =
  case o
  of boSingle, boNone: 1
  of boDouble: 2
  of boTriple: 3
  of boQuadruple: fail("quadruple bonds can't be written to a molfile", 0)
  of boAromatic: 4

func coord(m: Mol; i: int): Point3 =
  if m.coords.len == 0: (0.0, 0.0, 0.0) else: m.coords[i]

proc header(name, dim: string): string =
  ## Line 2 (CTfile): initials (2), program (8), date (10), dimension (2).
  name.replace("\n", " ") & "\n" & "  " & alignLeft("chem", 8) & spaces(10) & dim & "\n\n"

func fixed(x: float; width: int): string =
  ## `x` with 4 decimals, right-aligned in `width` columns.
  align(formatFloat(x, ffDecimal, 4), width)

func col(x: int; width = 3): string = align($x, width)

proc writeV2000(p: Prepared; name: string): string =
  let m = p.m
  result = header(name, m.dimension)
  result.add col(m.numAtoms) & col(m.numBonds) & "  0  0  0  0  0  0  0  0999 V2000\n"
  for i, a in m.atoms:
    let c = m.coord(i)
    let v = case writtenValence(m, i)
      of -1: 0
      of 0: 15
      else: writtenValence(m, i)
    result.add fixed(c.x, 10) & fixed(c.y, 10) & fixed(c.z, 10) & " " &
               alignLeft(symbol(a.element), 3) & " 0  0" & col(mdlParity(m, i)) &
               "  0  0" & col(v) & "  0  0  0" & col(a.atomMap) & "  0  0\n"
  const stereoCode = [bcNone: 0, bcWedge: 1, bcHash: 6, bcEither: 3]
  for bi, bd in m.bonds:
    let (a, b) = p.ends[bi]
    result.add col(a + 1) & col(b + 1) & col(bondTypeCode(bd.order)) &
               col(stereoCode[p.cfg[bi]]) & "\n"
  template nonzero(field): seq[(int, int)] =
    toSeq(0..<m.numAtoms).filterIt(m.atoms[it].field != 0).mapIt((it, m.atoms[it].field))
  for (key, values) in [("CHG", nonzero(charge)), ("ISO", nonzero(isotope))]:
    for start in countup(0, values.high, 8):
      let chunk = values[start .. min(start + 7, values.high)]
      result.add "M  " & key & col(chunk.len)
      for (i, v) in chunk: result.add " " & col(i + 1) & " " & col(v)
      result.add "\n"
  result.add "M  END\n"

proc writeV3000(p: Prepared; name: string): string =
  let m = p.m
  result = header(name, m.dimension)
  result.add "  0  0  0     0  0            999 V3000\n"
  result.add "M  V30 BEGIN CTAB\n"
  result.add "M  V30 COUNTS " & $m.numAtoms & " " & $m.numBonds & " 0 0 0\n"
  result.add "M  V30 BEGIN ATOM\n"
  for i, a in m.atoms:
    let c = m.coord(i)
    result.add "M  V30 " & $(i + 1) & " " & symbol(a.element) & " " &
               formatFloat(c.x, ffDecimal, 4) & " " & formatFloat(c.y, ffDecimal, 4) & " " &
               formatFloat(c.z, ffDecimal, 4) & " " & $a.atomMap
    if a.charge != 0: result.add " CHG=" & $a.charge
    if a.isotope != 0: result.add " MASS=" & $a.isotope
    let parity = mdlParity(m, i)
    if parity != 0: result.add " CFG=" & $parity
    let v = writtenValence(m, i)
    if v >= 0: result.add " VAL=" & $(if v == 0: -1 else: v)
    result.add "\n"
  result.add "M  V30 END ATOM\n"
  result.add "M  V30 BEGIN BOND\n"
  const cfgCode = [bcNone: 0, bcWedge: 1, bcHash: 3, bcEither: 2]
  for bi, bd in m.bonds:
    let (a, b) = p.ends[bi]
    result.add "M  V30 " & $(bi + 1) & " " & $bondTypeCode(bd.order) & " " &
               $(a + 1) & " " & $(b + 1)
    if p.cfg[bi] != bcNone: result.add " CFG=" & $cfgCode[p.cfg[bi]]
    result.add "\n"
  result.add "M  V30 END BOND\n"
  result.add "M  V30 END CTAB\n"
  result.add "M  END\n"

proc toMolBlock*(m: Mol; name = ""; v3000 = false; kekule = true): string =
  ## Writes `m` as a molfile: V2000, or V3000 when asked or when `m` has
  ## more than 999 atoms or bonds. With `kekule` (default) aromatic bonds
  ## are written as single and double bonds (as RDKit does), if possible.
  let p = prepare(m, kekule)
  if v3000 or m.numAtoms > 999 or m.numBonds > 999: writeV3000(p, name)
  else: writeV2000(p, name)
