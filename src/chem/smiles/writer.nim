## SMILES writer.
##
## Atoms are written in depth-first order, starting each connected component
## at its lowest-ranked atom and visiting neighbours in rank order. With
## canonical ranks (see `canon.nim`) the output is a canonical SMILES.

import std/[algorithm, sequtils, strutils]
import ../[mol, elements, valence, canon, rings, aromaticity]
import ./chirality

type Writer = object
  ## Traversal state. The molecule itself is passed separately (read-only)
  ## so it is never copied.
  sortedAdj: seq[seq[(int, int)]] ## adjacency sorted by neighbour rank
  inRing: seq[bool] ## per bond
  visited: seq[bool]
  bondUsed: seq[bool]
  parent: seq[int]
  children: seq[seq[(int, int)]] ## (child atom, bond) in output order
  openings: seq[seq[int]] ## ring bonds opened at this atom
  closings: seq[seq[int]] ## ring bonds closed at this atom
  ringDigit: seq[int] ## per bond
  digitInUse: array[1 .. 99, bool]
  pos: seq[int] ## per atom: output position
  writtenFrom: seq[int] ## per bond: atom it is written from
  mark: seq[string] ## per bond: "", "/" or "\\" (as written)
  treeBond: seq[bool] ## per bond: spanning-tree bond (not a ring closure)
  output: string

proc buildTree(w: var Writer, root: int) =
  ## Depth-first spanning tree from `root`, visiting neighbours in rank
  ## order. Iterative (explicit stack) so long chains can't overflow the
  ## call stack.
  w.visited[root] = true
  var stack = @[(root, 0)] # (atom, next index into sortedAdj)
  while stack.len > 0:
    let (a, idx) = stack[^1]
    if idx == w.sortedAdj[a].len:
      discard stack.pop
      continue
    stack[^1][1] = idx + 1
    let (nb, bi) = w.sortedAdj[a][idx]
    if w.bondUsed[bi]:
      continue
    w.bondUsed[bi] = true
    if w.visited[nb]:
      # back edge to an ancestor: ring closure, opened at nb, closed here
      w.openings[nb].add bi
      w.closings[a].add bi
    else:
      w.visited[nb] = true
      w.parent[nb] = a
      w.children[a].add (nb, bi)
      w.treeBond[bi] = true
      stack.add (nb, 0)

func atomToken(m: Mol, i: int, chiral: ChiralTag): string =
  let a = m.atoms[i]
  let hcount = totalHCount(m, i)
  let defaultH = defaultImplicitHCount(m, i)
  var sym = symbol(a.element)
  if a.aromatic:
    sym = sym.toLowerAscii
  let canShort =
    a.charge == 0 and a.isotope == 0 and chiral == ctNone and a.atomMap == 0 and (
      a.element == 0 or
      (if a.aromatic: inOrganicAromatic(a.element)
      else: inOrganicSubset(a.element))
    ) and (if a.element == 0: hcount == 0 else: hcount == defaultH)
  if canShort:
    return sym
  result = "["
  if a.isotope > 0:
    result &= $a.isotope
  result &= sym
  case chiral
  of ctCCW:
    result &= "@"
  of ctCW:
    result &= "@@"
  of ctNone:
    discard
  if hcount == 1:
    result &= "H"
  elif hcount > 1:
    result &= "H" & $hcount
  result &= chargeSuffix(a.charge)
  if a.atomMap != 0:
    result &= ":" & $a.atomMap
  result &= "]"

func bondToken(m: Mol, bi: int, inRing: bool, mark: string): string =
  ## `inRing` matters for aromatic bonds: an unwritten bond between aromatic
  ## atoms is only read back as aromatic inside a ring. `mark` is a '/' or
  ## '\\' double-bond stereo mark (only ever on single or aromatic bonds),
  ## already in written direction.
  if mark.len > 0:
    return mark
  let bd = m.bonds[bi]
  let bothArom = m.atoms[bd.a].aromatic and m.atoms[bd.b].aromatic
  case bd.order
  of boSingle:
    if bothArom: "-" else: ""
  of boDouble:
    "="
  of boTriple:
    "#"
  of boQuadruple:
    "$"
  of boAromatic:
    if bothArom and inRing: "" else: ":"
  of boNone:
    raise newException(
      ChemError,
      "cannot write bond " & $bi & " with unspecified order (boNone) as SMILES",
    )

proc emitBond(w: Writer, m: Mol, bi: int): string =
  bondToken(m, bi, w.inRing[bi], w.mark[bi])

proc layout(w: var Writer, root: int): seq[int] =
  ## Records output positions (pre-order, as `writeComponent` will write)
  ## and the direction each bond is written in. Returns the atoms in order.
  var stack = @[root]
  while stack.len > 0:
    let a = stack.pop
    w.pos[a] = result.len
    result.add a
    for bi in w.openings[a]:
      w.writtenFrom[bi] = a # ring bond: from opener
    for k in countdown(w.children[a].high, 0):
      let (c, bi) = w.children[a][k]
      w.writtenFrom[bi] = a
      stack.add c

proc assignStereoMarks(w: var Writer, m: Mol, atoms: seq[int]) =
  ## Chooses '/' '\\' marks for the stereo double bonds of one component.
  ##
  ## 1. Each end of each stereo double bond gets a carrier bond (a single or
  ##    aromatic bond to a substituent), reusing one already chosen for a
  ##    conjugated neighbour, else the first in output order (tree bonds
  ##    before ring closures).
  ## 2. Every double bond then constrains the marks on all carriers at its
  ##    ends (substituents on one end point opposite ways; across the bond
  ##    they match for cis, differ for trans). These are parity constraints,
  ##    solved by 2-colouring with the first-written carrier of each group
  ##    set to '/', so equivalent inputs give identical output.
  var dbs: seq[(int, int)] # (position of first-written end, bond)
  for a in atoms:
    for (nb, bi) in m.adj[a]:
      if m.bonds[bi].stereo != bsNone and w.pos[a] < w.pos[nb]:
        dbs.add (w.pos[a], bi)
  if dbs.len == 0:
    return
  dbs.sort()

  proc carrierBonds(e, other: int): seq[(int, int)] =
    for (nb, bi) in m.adj[e]:
      if nb != other and m.bonds[bi].order in {boSingle, boAromatic}:
        result.add (nb, bi)

  # 1. carriers
  var chosen = newSeq[bool](m.numBonds)
  for (_, bi) in dbs:
    let bd = m.bonds[bi]
    for (e, other) in [(bd.a, bd.b), (bd.b, bd.a)]:
      let cs = carrierBonds(e, other)
      if cs.anyIt(chosen[it[1]]):
        continue
      var best: (bool, int, int) # (ring closure, position, bond)
      best[2] = -1
      for (nb, cb) in cs:
        let key = (not w.treeBond[cb], w.pos[nb], cb)
        if best[2] < 0 or key < best:
          best = key
      if best[2] < 0:
        raise newException(
          ChemError, "no bond can carry the '/' '\\' mark for double bond " & $bi
        )
      chosen[best[2]] = true

  # 2. parity constraints: value[c] is 0 for '/', 1 for '\\' as written;
  # a mark read outward from atom e is value xor (written from the other end)
  var edges = newSeq[seq[(int, int)]](m.numBonds) # (other bond, parity)
  proc addConstraint(c1, c2, parity: int) =
    edges[c1].add (c2, parity)
    edges[c2].add (c1, parity)

  for (_, bi) in dbs:
    let bd = m.bonds[bi]
    var ends: array[2, seq[(int, int, int)]] # (substituent, bond, flipped)
    for k, (e, other) in [(bd.a, bd.b), (bd.b, bd.a)]:
      for (nb, cb) in carrierBonds(e, other):
        if chosen[cb]:
          ends[k].add (nb, cb, ord(w.writtenFrom[cb] != e))
    for k in 0 .. 1: # two substituents on one end point opposite ways
      if ends[k].len == 2:
        let (_, c1, f1) = ends[k][0]
        let (_, c2, f2) = ends[k][1]
        addConstraint(c1, c2, 1 xor f1 xor f2)
    for (n1, c1, f1) in ends[0]:
      for (n2, c2, f2) in ends[1]:
        let cis =
          (bd.stereo == bsCis) xor (n1 != bd.stereoAtoms[0]) xor
          (n2 != bd.stereoAtoms[1])
        addConstraint(c1, c2, ord(not cis) xor f1 xor f2)

  # 3. 2-colouring, groups rooted at their first-written carrier
  var order: seq[(int, int, int)] # (position written from, other end, bond)
  for cb in 0 ..< m.numBonds:
    if chosen[cb]:
      let bd = m.bonds[cb]
      let src = w.writtenFrom[cb]
      order.add (w.pos[src], w.pos[bd.otherAtom(src)], cb)
  order.sort()
  var value = newSeqWith(m.numBonds, -1)
  for (_, _, root) in order:
    if value[root] >= 0:
      continue
    value[root] = 0
    var queue = @[root]
    while queue.len > 0:
      let c = queue.pop
      for (d, parity) in edges[c]:
        let want = value[c] xor parity
        if value[d] < 0:
          value[d] = want
          queue.add d
        elif value[d] != want:
          raise newException(
            ChemError,
            "cannot write consistent '/' '\\' marks for the double bonds around bond " &
              $c,
          )
  for cb in 0 ..< m.numBonds:
    if value[cb] >= 0:
      w.mark[cb] = if value[cb] == 0: "/" else: "\\"

func digitToken(d: int): string =
  if d < 10:
    $d
  else:
    "%" & $d

proc writeAtomToken(w: var Writer, m: Mol, a: int) =
  ## Writes atom `a` and its ring-closure digits (not its branches).
  # re-express the chirality tag for the neighbour order as written
  var chiral = m.atoms[a].chiral
  if chiral != ctNone:
    let reference = chiralReference(m, a)
    var written: seq[int]
    if w.parent[a] >= 0:
      written.add w.parent[a]
    if reference.len > 0 and reference[0] == -1:
      written.add -1
    for bi in w.closings[a]:
      written.add m.bonds[bi].otherAtom(a)
    for bi in w.openings[a]:
      written.add m.bonds[bi].otherAtom(a)
    for (c, _) in w.children[a]:
      written.add c
    chiral = reorderChirality(
      m, a, chiral, reference, written, w.closings[a].len + w.openings[a].len
    )
  w.output &= atomToken(m, a, chiral)
  var toFree: seq[int]
  for bi in w.closings[a]:
    let d = w.ringDigit[bi]
    w.output &= digitToken(d)
    toFree.add d
  for bi in w.openings[a]:
    var d = 1
    while d <= 99 and w.digitInUse[d]:
      inc d
    if d > 99:
      raise newException(ChemError, "too many open rings for SMILES output")
    w.digitInUse[d] = true
    w.ringDigit[bi] = d
    w.output &= w.emitBond(m, bi) & digitToken(d)
  for d in toFree:
    w.digitInUse[d] = false

proc writeComponent(w: var Writer, m: Mol, root: int) =
  ## Writes the spanning tree from `root`: all but the last child of each
  ## atom become parenthesised branches. Iterative, like `buildTree`.
  w.writeAtomToken(m, root)
  var stack = @[(root, 0)] # (atom, next child index)
  while stack.len > 0:
    let (a, k) = stack[^1]
    let last = w.children[a].high
    if k > 0 and k - 1 < last:
      w.output &= ")" # close previous branch
    if k > last:
      discard stack.pop
      continue
    stack[^1][1] = k + 1
    let (c, bi) = w.children[a][k]
    if k < last:
      w.output &= "("
    w.output &= w.emitBond(m, bi)
    w.writeAtomToken(m, c)
    stack.add (c, 0)

proc toSmilesWithRanks*(m: Mol, ranks: seq[int]): string =
  ## Writes SMILES using the given atom ranks to choose the traversal order.
  let n = m.numAtoms
  var w: Writer
  w.sortedAdj = m.adj
  for nbrs in w.sortedAdj.mitems:
    nbrs.sort(
      proc(x, y: (int, int)): int =
        cmp(ranks[x[0]], ranks[y[0]])
    )
  w.inRing = ringBonds(m)
  w.visited = newSeq[bool](n)
  w.bondUsed = newSeq[bool](m.numBonds)
  w.parent = newSeqWith(n, -1)
  w.children = newSeq[seq[(int, int)]](n)
  w.openings = newSeq[seq[int]](n)
  w.closings = newSeq[seq[int]](n)
  w.ringDigit = newSeq[int](m.numBonds)
  w.pos = newSeq[int](n)
  w.writtenFrom = newSeqWith(m.numBonds, -1)
  w.mark = newSeq[string](m.numBonds)
  w.treeBond = newSeq[bool](m.numBonds)
  # components ordered by their lowest-ranked atom
  let comp = connectedComponents(m)
  var starts = newSeqWith(componentCount(comp), -1)
  for i in 0 ..< n:
    let c = comp[i]
    if starts[c] < 0 or ranks[i] < ranks[starts[c]]:
      starts[c] = i
  starts.sort(
    proc(x, y: int): int =
      cmp(ranks[x], ranks[y])
  )
  var parts: seq[string]
  for s in starts:
    w.buildTree(s)
    let atoms = w.layout(s)
    w.assignStereoMarks(m, atoms)
    w.output = ""
    w.writeComponent(m, s)
    parts.add w.output
  parts.join(".")

proc toSmiles*(m: Mol, canonical = true, kekule = false): string

proc canonicalByComponent(m: Mol, comp: seq[int], kekule: bool): string =
  ## Canonical SMILES of each component on its own, sorted. Canonicalising
  ## the whole molecule at once can order identical-looking components
  ## (e.g. cis and trans isomers in one mixture) differently depending on
  ## input order.
  var members = newSeq[seq[int]](componentCount(comp))
  for i, c in comp:
    members[c].add i
  var parts = members.mapIt(m.subMol(it).toSmiles(canonical = true, kekule = kekule))
  parts.sort()
  parts.join(".")

proc toSmiles*(m: Mol, canonical = true, kekule = false): string =
  ## Writes a SMILES string. With `canonical` (default) the result is
  ## independent of the input atom order. With `kekule`, aromatic systems
  ## are written as alternating single/double bonds (uppercase atoms).
  if canonical:
    let comp = connectedComponents(m)
    if componentCount(comp) > 1:
      return canonicalByComponent(m, comp, kekule)
  if kekule:
    # kekulize in canonical atom order so the double-bond placement is
    # canonical too
    var k =
      if canonical:
        let ranks = canonicalRanks(m)
        var order = newSeq[int](m.numAtoms)
        for i, r in ranks:
          order[r] = i
        m.renumberAtoms(order)
      else:
        m
    kekulize(k)
    toSmilesWithRanks(k, toSeq(0 ..< k.numAtoms))
  elif canonical:
    toSmilesWithRanks(m, canonicalRanks(m))
  else:
    toSmilesWithRanks(m, toSeq(0 ..< m.numAtoms))
