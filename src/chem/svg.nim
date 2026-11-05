## Drawing molecules as SVG: bonds (Kekulé form, with double bonds inside
## rings), wedge and hash bonds for tetrahedral stereo, crossed double
## bonds where cis/trans is unknown, atom labels (heteroatoms with their
## hydrogens, charges, isotopes), element colours, highlighting, and grids
## of molecules with legends. Molecules without 2D coordinates are laid out
## first (see `depict`); 3D coordinates are laid out afresh too.

import std/[math, strutils, sequtils, unicode]
import ./[mol, elements, valence, rings, aromaticity, geometry, depict, molfile]

type
  DrawOptions* = object
    width* = 0 ## image size in px; 0 fits the drawing
    height* = 0
    bondLength* = 30.0 ## px (shrunk if needed to fit a given size)
    fontSize* = 0.0 ## px; 0 is 0.45 bond lengths
    lineWidth* = 0.0 ## px; 0 is 1/15 bond length
    padding* = 10.0 ## px around the drawing
    background* = "#ffffff" ## "" for transparent
    colorAtoms* = true ## element colours for labels and bonds
    highlightColor* = "#ff8080"

  Pt = tuple[x, y: float]

  Drawing = object ## A molecule ready to draw: positions in px, y downwards.
    m: Mol
    pos: seq[Pt]
    labelled: seq[bool]
    wedges: seq[(int, BondWedge)]
    either: seq[bool]
    ringCentre: seq[(bool, Pt)] ## per bond: inside which ring (centre)
    bl, font, lw: float

func `+`(a, b: Pt): Pt =
  (a.x + b.x, a.y + b.y)
func `-`(a, b: Pt): Pt =
  (a.x - b.x, a.y - b.y)
func `*`(a: Pt, k: float): Pt =
  (a.x * k, a.y * k)
func len(a: Pt): float =
  hypot(a.x, a.y)
func unit(a: Pt): Pt =
  let l = len(a)
  if l == 0:
    (0.0, 0.0)
  else:
    a * (1 / l)
func normal(a: Pt): Pt =
  unit((-a.y, a.x))

func f(x: float): string =
  formatFloat(x, ffDecimal, 1)

func xmlEscape(s: string): string =
  for c in s:
    case c
    of '&':
      result.add "&amp;"
    of '<':
      result.add "&lt;"
    of '>':
      result.add "&gt;"
    of '"':
      result.add "&quot;"
    else:
      result.add c

func elementColor(z: int): string =
  ## colours as RDKit's default palette
  case z
  of 7: "#0000ff"
  of 8: "#ff0000"
  of 9: "#33cccc"
  of 15: "#ff8000"
  of 16: "#cccc00"
  of 17: "#00cc00"
  of 35: "#994d1a"
  of 53: "#a01ef0"
  of 1, 0: "#7f7f7f"
  else: "#000000"

# ---------------------------------------------------------------------------
# Preparing

func needsLabel(m: Mol, i: int): bool =
  let a = m.atoms[i]
  a.element != 6 or a.charge != 0 or a.isotope != 0 or a.atomMap != 0 or m.degree(i) == 0

proc prepareDrawing(m: Mol, opts: DrawOptions, bl: float): Drawing =
  var k = m
  if k.coords.len == 0 or k.is3D:
    compute2DCoords(k)
  try:
    kekulize(k)
  except ChemError:
    discard # aromatic bonds are drawn dashed
  result.m = k
  result.bl = bl
  result.font =
    if opts.fontSize > 0:
      opts.fontSize
    else:
      0.45 * bl
  result.lw =
    if opts.lineWidth > 0:
      opts.lineWidth
    else:
      bl / 15
  # scale to `bl` px per (mean) bond, y downwards
  var mean = 0.0
  for bd in k.bonds:
    mean +=
      hypot(k.coords[bd.a].x - k.coords[bd.b].x, k.coords[bd.a].y - k.coords[bd.b].y)
  mean =
    if k.numBonds > 0 and mean > 0:
      mean / float(k.numBonds)
    else:
      bondLength
  let s = bl / mean
  result.pos = k.coords.mapIt((it.x * s, -it.y * s))
  result.labelled = toSeq(0 ..< k.numAtoms).mapIt(needsLabel(k, it))
  try:
    result.wedges = chooseWedges(k)
  except ChemError:
    result.wedges = newSeqWith(k.numBonds, (-1, bwNone))
  result.either = unspecifiedStereoBonds(k)
  # which ring each ring double bond is drawn inside: the smallest, and
  # among equals one whose bonds are all double/aromatic-ish
  result.ringCentre = newSeq[(bool, Pt)](k.numBonds)
  var ringSize = newSeqWith(k.numBonds, high(int))
  for ring in sssr(k):
    var c: Pt
    for x in ring:
      c = c + result.pos[x]
    c = c * (1 / float(ring.len))
    for bi in k.ringBondIndices(ring):
      if ring.len < ringSize[bi]:
        ringSize[bi] = ring.len
        result.ringCentre[bi] = (true, c)

# ---------------------------------------------------------------------------
# Labels

type Label = object ## Text pieces around an atom, positioned relative to it.
  parts: seq[(float, string, string)] ## (x offset, anchor, svg content)
  halfWidth: float ## how far bonds stop from the atom

proc atomLabel(d: Drawing, i: int): Label =
  let a = d.m.atoms[i]
  let font = d.font
  let cw = 0.6 * font # estimated character width
  let sym =
    if a.element == 0:
      "*"
    else:
      symbol(a.element)
  let h = totalHCount(d.m, i)
  # H on the left when the bonds leave to the right
  var dx = 0.0
  for nb in d.m.neighbors(i):
    dx += d.pos[nb].x - d.pos[i].x
  let hLeft = dx > 1e-6 and d.m.degree(i) > 0
  let symW = cw * float(sym.len)
  result.halfWidth = max(symW, cw) * 0.55 + 0.1 * font
  result.parts.add (0.0, "middle", xmlEscape(sym))
  let small = "font-size=\"" & f(0.7 * font) & "\""
  var hText = ""
  if h > 0:
    hText = "H"
    if h > 1:
      hText &= "<tspan " & small & " dy=\"" & f(0.25 * font) & "\">" & $h & "</tspan>"
  var charge = ""
  if a.charge != 0:
    let n = abs(a.charge)
    charge = (if n > 1: $n else: "") & (if a.charge > 0: "+" else: "−")
  let chargeSpan =
    if charge.len > 0:
      "<tspan " & small & " dy=\"" & f(-0.4 * font) & "\">" & charge & "</tspan>"
    else:
      ""
  let hW =
    if h > 1:
      cw * 1.7
    else:
      cw
  # the isotope goes just before the symbol, H on the left before that
  let isoW =
    if a.isotope > 0:
      0.7 * cw * float(len($a.isotope))
    else:
      0.0
  if a.isotope > 0:
    result.parts.add (
      -symW / 2,
      "end",
      "<tspan " & small & " dy=\"" & f(-0.4 * font) & "\">" & $a.isotope & "</tspan>",
    )
  if hLeft and hText.len > 0:
    result.parts.add (-symW / 2 - isoW, "end", hText)
    if charge.len > 0:
      result.parts.add (symW / 2, "start", chargeSpan)
  else:
    var tail = hText & chargeSpan
    if tail.len > 0:
      result.parts.add (symW / 2, "start", tail)
  if a.atomMap > 0:
    result.parts.add (
      symW / 2 + (if not hLeft and h > 0: hW else: 0.0),
      "start",
      "<tspan " & small & " dy=\"" & f(0.35 * font) & "\">:" & $a.atomMap & "</tspan>",
    )

# ---------------------------------------------------------------------------
# Drawing

proc color(d: Drawing, i: int, opts: DrawOptions): string =
  if opts.colorAtoms:
    elementColor(d.m.atoms[i].element)
  else:
    "#000000"

proc line(o: var string, a, b: Pt, col: string, lw: float, dash = false) =
  o.add "<path d=\"M " & f(a.x) & "," & f(a.y) & " L " & f(b.x) & "," & f(b.y) &
    "\" stroke=\"" & col & "\" stroke-width=\"" & f(lw) & "\" stroke-linecap=\"round\"" &
    (if dash: " stroke-dasharray=\"" & f(lw * 2) & "," & f(lw * 2) & "\"" else: "") &
    "/>\n"

proc halfLines(
    o: var string, d: Drawing, a, b: Pt, ia, ib: int, opts: DrawOptions, dash = false
) =
  ## A line from `a` to `b` coloured half by atom `ia`'s colour, half by `ib`'s.
  let ca = d.color(ia, opts)
  let cb = d.color(ib, opts)
  if ca == cb:
    line(o, a, b, ca, d.lw, dash)
  else:
    let mid = (a + b) * 0.5
    line(o, a, mid, ca, d.lw, dash)
    line(o, mid, b, cb, d.lw, dash)

proc drawBond(o: var string, d: Drawing, bi: int, opts: DrawOptions) =
  let bd = d.m.bonds[bi]
  var (i, j) = (bd.a, bd.b)
  let (narrow, wedge) = d.wedges[bi]
  if wedge != bwNone and narrow == j:
    swap(i, j)
  var a = d.pos[i]
  var b = d.pos[j]
  let u = unit(b - a)
  let labels = [d.labelled[i], d.labelled[j]]
  # stop short of labels
  let ri =
    if labels[0]:
      atomLabel(d, i).halfWidth
    else:
      0.0
  let rj =
    if labels[1]:
      atomLabel(d, j).halfWidth
    else:
      0.0
  if len(b - a) > ri + rj:
    a = a + u * ri
    b = b - u * rj
  let nrm = normal(b - a)
  let gap = 0.18 * d.bl
  case bd.order
  of boSingle, boNone, boQuadruple, boAromatic:
    if wedge == bwWedge:
      let w = 0.12 * d.bl
      o.add "<path d=\"M " & f(a.x) & "," & f(a.y) & " L " & f((b + nrm * w).x) & "," &
        f((b + nrm * w).y) & " L " & f((b - nrm * w).x) & "," & f((b - nrm * w).y) &
        " Z\" fill=\"" & d.color(i, opts) & "\" stroke=\"" & d.color(i, opts) &
        "\" stroke-width=\"" & f(d.lw * 0.5) & "\"/>\n"
    elif wedge == bwHash:
      let w = 0.12 * d.bl
      const n = 6
      for k in 1 .. n:
        let t = float(k) / float(n)
        let p = a + (b - a) * t
        line(o, p + nrm * (w * t), p - nrm * (w * t), d.color(i, opts), d.lw * 0.8)
    else:
      halfLines(o, d, a, b, i, j, opts)
    if bd.order == boAromatic or bd.order == boQuadruple:
      # a dashed (aromatic) or extra line inside the ring / to one side
      let (inRing, c) = d.ringCentre[bi]
      var side = nrm
      if inRing and (c - a).x * nrm.x + (c - a).y * nrm.y < 0:
        side = nrm * -1.0
      let shrink = (b - a) * 0.15
      halfLines(
        o,
        d,
        a + side * gap + shrink,
        b + side * gap - shrink,
        i,
        j,
        opts,
        dash = bd.order == boAromatic,
      )
  of boDouble:
    if d.either[bi]:
      let h = nrm * (gap / 2)
      halfLines(o, d, a + h, b - h, i, j, opts)
      halfLines(o, d, a - h, b + h, i, j, opts)
      return
    let (inRing, c) = d.ringCentre[bi]
    let terminal = d.m.degree(i) == 1 or d.m.degree(j) == 1
    if inRing or not terminal:
      # the main line, and a shorter one inside the ring (or on the side
      # with more substituents)
      var side = nrm
      if inRing:
        if (c - a).x * nrm.x + (c - a).y * nrm.y < 0:
          side = nrm * -1.0
      else:
        var score = 0.0
        for x in [i, j]:
          for nb in d.m.neighbors(x):
            if nb == i or nb == j:
              continue
            score += (
              if (d.pos[nb] - a).x * nrm.x + (d.pos[nb] - a).y * nrm.y > 0: 1.0
              else: -1.0
            )
        if score < 0:
          side = nrm * -1.0
      halfLines(o, d, a, b, i, j, opts)
      let sa = if labels[0]: 0.0 else: 0.15
      let sb = if labels[1]: 0.0 else: 0.15
      let ab = b - a
      halfLines(o, d, a + side * gap + ab * sa, b + side * gap - ab * sb, i, j, opts)
    else:
      let h = nrm * (gap / 2)
      halfLines(o, d, a + h, b + h, i, j, opts)
      halfLines(o, d, a - h, b - h, i, j, opts)
  of boTriple:
    halfLines(o, d, a, b, i, j, opts)
    halfLines(o, d, a + nrm * gap, b + nrm * gap, i, j, opts)
    halfLines(o, d, a - nrm * gap, b - nrm * gap, i, j, opts)

proc bounds(d: Drawing): (Pt, Pt) =
  ## extent of the atoms, with room for labels
  var lo: Pt = (Inf, Inf)
  var hi: Pt = (-Inf, -Inf)
  for i, p in d.pos:
    let r =
      if d.labelled[i]:
        d.font * 1.6
      else:
        0.0
    lo = (min(lo.x, p.x - r), min(lo.y, p.y - d.font * 0.7))
    hi = (max(hi.x, p.x + r), max(hi.y, p.y + d.font * 0.7))
  if d.pos.len == 0:
    return ((0.0, 0.0), (0.0, 0.0))
  (lo, hi)

proc drawBody(
    m: Mol,
    opts: DrawOptions,
    w, h: float,
    highlightAtoms, highlightBonds: openArray[int],
): (string, float, float) =
  ## The drawing's elements, fitted into `w` x `h` px (0: natural size),
  ## and the size used.
  var d = prepareDrawing(m, opts, opts.bondLength)
  var (lo, hi) = bounds(d)
  var (width, height) = (hi.x - lo.x + 2 * opts.padding, hi.y - lo.y + 2 * opts.padding)
  if w > 0 and h > 0:
    # shrink to fit (never enlarge beyond the bond length)
    let k = min(
      1.0,
      min(
        (w - 2 * opts.padding) / max(hi.x - lo.x, 1e-9),
        (h - 2 * opts.padding) / max(hi.y - lo.y, 1e-9),
      ),
    )
    if k < 1:
      d = prepareDrawing(m, opts, opts.bondLength * k)
      (lo, hi) = bounds(d)
    (width, height) = (w, h)
  let shift: Pt =
    ((width - (hi.x - lo.x)) / 2 - lo.x, (height - (hi.y - lo.y)) / 2 - lo.y)
  for p in d.pos.mitems:
    p = p + shift
  for c in d.ringCentre.mitems:
    c[1] = c[1] + shift
  var o = ""
  # highlights first, underneath
  for bi in highlightBonds:
    if bi in 0 ..< d.m.numBonds:
      line(
        o,
        d.pos[d.m.bonds[bi].a],
        d.pos[d.m.bonds[bi].b],
        opts.highlightColor,
        0.35 * d.bl,
      )
  for i in highlightAtoms:
    if i in 0 ..< d.m.numAtoms:
      o.add "<circle cx=\"" & f(d.pos[i].x) & "\" cy=\"" & f(d.pos[i].y) & "\" r=\"" &
        f(0.3 * d.bl) & "\" fill=\"" & opts.highlightColor & "\"/>\n"
  for bi in 0 ..< d.m.numBonds:
    drawBond(o, d, bi, opts)
  for i in 0 ..< d.m.numAtoms:
    if not d.labelled[i]:
      continue
    let lab = atomLabel(d, i)
    let p = d.pos[i]
    for (dx, anchor, content) in lab.parts:
      o.add "<text x=\"" & f(p.x + dx) & "\" y=\"" & f(p.y + d.font * 0.35) &
        "\" text-anchor=\"" & anchor & "\" fill=\"" & d.color(i, opts) & "\">" & content &
        "</text>\n"
  (o, width, height)

func svgHeader(w, h: float, opts: DrawOptions, font: float): string =
  result =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" &
    "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"" & f(w) &
    "\" height=\"" & f(h) & "\" viewBox=\"0 0 " & f(w) & " " & f(h) & "\" " &
    "font-family=\"sans-serif\" font-size=\"" & f(font) & "\">\n"
  if opts.background.len > 0:
    result.add "<rect width=\"100%\" height=\"100%\" fill=\"" & opts.background &
      "\"/>\n"

proc toSvg*(
    m: Mol,
    opts = DrawOptions(),
    highlightAtoms: openArray[int] = [],
    highlightBonds: openArray[int] = [],
): string =
  ## An SVG drawing of `m`. Highlighted atoms and bonds (e.g. a substructure
  ## match, with `highlightBonds` the bonds between its atoms) get a
  ## coloured background.
  let font =
    if opts.fontSize > 0:
      opts.fontSize
    else:
      0.45 * opts.bondLength
  let (body, w, h) = drawBody(
    m, opts, float(opts.width), float(opts.height), highlightAtoms, highlightBonds
  )
  svgHeader(w, h, opts, font) & body & "</svg>\n"

func bondsBetween*(m: Mol, atoms: openArray[int]): seq[int] =
  ## The bonds joining two of `atoms` (to highlight with them).
  for bi, bd in m.bonds:
    if bd.a in atoms and bd.b in atoms:
      result.add bi

proc toSvgGrid*(
    mols: openArray[Mol],
    legends: openArray[string] = [],
    perRow = 4,
    cellWidth = 250,
    cellHeight = 200,
    opts = DrawOptions(),
    highlights: openArray[seq[int]] = [],
): string =
  ## Several molecules in a grid of `perRow` columns, each in a
  ## `cellWidth` x `cellHeight` px cell (drawings shrunk to fit), with an
  ## optional legend under each, and optionally atoms to highlight in each
  ## (with the bonds between them).
  let rows = max(1, (mols.len + perRow - 1) div perRow)
  let cols = max(1, min(perRow, mols.len))
  let font =
    if opts.fontSize > 0:
      opts.fontSize
    else:
      0.45 * opts.bondLength
  let legendH =
    if legends.len > 0:
      font * 1.5
    else:
      0.0
  let (w, h) = (float(cols * cellWidth), float(rows * cellHeight))
  result = svgHeader(w, h, opts, font)
  for k, m in mols:
    let (x, y) = (float((k mod perRow) * cellWidth), float((k div perRow) * cellHeight))
    let hl =
      if k < highlights.len:
        highlights[k]
      else:
        @[]
    let (body, _, _) = drawBody(
      m, opts, float(cellWidth), float(cellHeight) - legendH, hl, m.bondsBetween(hl)
    )
    result.add "<g transform=\"translate(" & f(x) & "," & f(y) & ")\">\n" & body
    if k < legends.len:
      # (shortened to fit the cell, by an estimated character width)
      let fits = max(1, int(float(cellWidth) / (0.6 * font)) - 1)
      let legend =
        if legends[k].runeLen <= fits:
          legends[k]
        else:
          legends[k].runeSubStr(0, fits - 1) & "\u2026"
      result.add "<text x=\"" & f(float(cellWidth) / 2) & "\" y=\"" &
        f(float(cellHeight) - legendH * 0.35) & "\" text-anchor=\"middle\">" &
        xmlEscape(legend) & "</text>\n"
    result.add "</g>\n"
  result.add "</svg>\n"
