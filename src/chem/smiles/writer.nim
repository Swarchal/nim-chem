## SMILES writer.
##
## Atoms are written in depth-first order, starting each connected component
## at its lowest-ranked atom and visiting neighbours in rank order. With
## canonical ranks (see `canon.nim`) the output is a canonical SMILES.

import std/[algorithm, sequtils, strutils]
import ../[mol, elements, valence, canon, rings, aromaticity]

type
  Writer = object
    ## Traversal state. The molecule itself is passed separately (read-only)
    ## so it is never copied.
    sortedAdj: seq[seq[(int, int)]]  ## adjacency sorted by neighbour rank
    inRing: seq[bool]                ## per bond
    visited: seq[bool]
    bondUsed: seq[bool]
    parent: seq[int]
    children: seq[seq[(int, int)]]   ## (child atom, bond) in output order
    openings: seq[seq[int]]          ## ring bonds opened at this atom
    closings: seq[seq[int]]          ## ring bonds closed at this atom
    ringDigit: seq[int]              ## per bond
    digitInUse: array[1..99, bool]
    dirGroup: seq[int]               ## per bond: directional-bond group, -1 if none
    groupFlip: seq[int]              ## per group: -1 undecided, 0 keep, 1 invert
    output: string

func directionalGroups(m: Mol): (seq[int], int) =
  ## Groups '/' '\' bonds that are coupled through double bonds; inverting
  ## every marker in a group leaves the double-bond geometry unchanged.
  var parent = toSeq(0..<m.numBonds)
  proc find(x: int): int =
    var x = x
    while parent[x] != x:
      parent[x] = parent[parent[x]]
      x = parent[x]
    x
  for bd in m.bonds:
    if bd.order != boDouble: continue
    var first = -1
    for endAtom in [bd.a, bd.b]:
      for (_, bi) in m.adj[endAtom]:
        if m.bonds[bi].stereo == bsNone: continue
        if first < 0: first = bi
        else: parent[find(bi)] = find(first)
  var ids = newSeqWith(m.numBonds, -1)
  var count = 0
  var rootId = newSeqWith(m.numBonds, -1)
  for bi, bd in m.bonds:
    if bd.stereo != bsNone:
      let r = find(bi)
      if rootId[r] < 0:
        rootId[r] = count
        inc count
      ids[bi] = rootId[r]
  (ids, count)

proc buildTree(w: var Writer; root: int) =
  ## Depth-first spanning tree from `root`, visiting neighbours in rank
  ## order. Iterative (explicit stack) so long chains can't overflow the
  ## call stack.
  w.visited[root] = true
  var stack = @[(root, 0)]   # (atom, next index into sortedAdj)
  while stack.len > 0:
    let (a, idx) = stack[^1]
    if idx == w.sortedAdj[a].len:
      discard stack.pop
      continue
    stack[^1][1] = idx + 1
    let (nb, bi) = w.sortedAdj[a][idx]
    if w.bondUsed[bi]: continue
    w.bondUsed[bi] = true
    if w.visited[nb]:
      # back edge to an ancestor: ring closure, opened at nb, closed here
      w.openings[nb].add bi
      w.closings[a].add bi
    else:
      w.visited[nb] = true
      w.parent[nb] = a
      w.children[a].add (nb, bi)
      stack.add (nb, 0)

func atomToken(m: Mol; i: int; chiral: ChiralTag): string =
  let a = m.atoms[i]
  let hcount = totalHCount(m, i)
  let defaultH = defaultImplicitHCount(m, i)
  var sym = symbol(a.element)
  if a.aromatic: sym = sym.toLowerAscii
  let canShort = a.charge == 0 and a.isotope == 0 and chiral == ctNone and
                 a.atomMap == 0 and
                 (a.element == 0 or
                  (if a.aromatic: inOrganicAromatic(a.element)
                   else: inOrganicSubset(a.element))) and
                 (if a.element == 0: hcount == 0
                  else: hcount == defaultH)
  if canShort: return sym
  result = "["
  if a.isotope > 0: result &= $a.isotope
  result &= sym
  case chiral
  of ctCCW: result &= "@"
  of ctCW: result &= "@@"
  of ctNone: discard
  if hcount == 1: result &= "H"
  elif hcount > 1: result &= "H" & $hcount
  result &= chargeSuffix(a.charge)
  if a.atomMap != 0: result &= ":" & $a.atomMap
  result &= "]"

func bondToken(m: Mol; bi, fromAtom: int; inRing: bool; invert = false): string =
  ## `inRing` matters for aromatic bonds: an unwritten bond between aromatic
  ## atoms is only read back as aromatic inside a ring.
  let bd = m.bonds[bi]
  let bothArom = m.atoms[bd.a].aromatic and m.atoms[bd.b].aromatic
  case bd.order
  of boSingle:
    if bd.stereo != bsNone:
      let up = ((bd.stereo == bsUp) == (bd.a == fromAtom)) xor invert
      if up: "/" else: "\\"
    elif bothArom: "-"
    else: ""
  of boDouble: "="
  of boTriple: "#"
  of boQuadruple: "$"
  of boAromatic: (if bothArom and inRing: "" else: ":")
  of boNone:
    raise newException(ChemError, "cannot write bond " & $bi &
                       " with unspecified order (boNone) as SMILES")

proc emitBond(w: var Writer; m: Mol; bi, fromAtom: int): string =
  ## Bond token; the first directional marker written in each group is
  ## normalised to '/' so equivalent inputs give identical output.
  let g = w.dirGroup[bi]
  let inRing = w.inRing[bi]
  if g < 0: return bondToken(m, bi, fromAtom, inRing)
  if w.groupFlip[g] < 0:
    w.groupFlip[g] = ord(bondToken(m, bi, fromAtom, inRing) == "\\")
  bondToken(m, bi, fromAtom, inRing, w.groupFlip[g] == 1)

func digitToken(d: int): string =
  if d < 10: $d else: "%" & $d

proc writeAtomToken(w: var Writer; m: Mol; a: int) =
  ## Writes atom `a` and its ring-closure digits (not its branches).
  # re-express the chirality tag for the neighbour order as written
  var chiral = m.atoms[a].chiral
  if chiral != ctNone:
    let reference = chiralReference(m, a)
    var written: seq[int]
    if w.parent[a] >= 0: written.add w.parent[a]
    if reference.len > 0 and reference[0] == -1: written.add -1
    for bi in w.closings[a]: written.add m.bonds[bi].otherAtom(a)
    for bi in w.openings[a]: written.add m.bonds[bi].otherAtom(a)
    for (c, _) in w.children[a]: written.add c
    if permutationParity(reference, written) == 1:
      chiral = flip(chiral)
  w.output &= atomToken(m, a, chiral)
  var toFree: seq[int]
  for bi in w.closings[a]:
    let d = w.ringDigit[bi]
    w.output &= digitToken(d)
    toFree.add d
  for bi in w.openings[a]:
    var d = 1
    while d <= 99 and w.digitInUse[d]: inc d
    if d > 99:
      raise newException(ChemError, "too many open rings for SMILES output")
    w.digitInUse[d] = true
    w.ringDigit[bi] = d
    w.output &= w.emitBond(m, bi, a) & digitToken(d)
  for d in toFree: w.digitInUse[d] = false

proc writeComponent(w: var Writer; m: Mol; root: int) =
  ## Writes the spanning tree from `root`: all but the last child of each
  ## atom become parenthesised branches. Iterative, like `buildTree`.
  w.writeAtomToken(m, root)
  var stack = @[(root, 0)]   # (atom, next child index)
  while stack.len > 0:
    let (a, k) = stack[^1]
    let last = w.children[a].high
    if k > 0 and k - 1 < last: w.output &= ")"   # close previous branch
    if k > last:
      discard stack.pop
      continue
    stack[^1][1] = k + 1
    let (c, bi) = w.children[a][k]
    if k < last: w.output &= "("
    w.output &= w.emitBond(m, bi, a)
    w.writeAtomToken(m, c)
    stack.add (c, 0)

proc toSmilesWithRanks*(m: Mol; ranks: seq[int]): string =
  ## Writes SMILES using the given atom ranks to choose the traversal order.
  let n = m.numAtoms
  var w: Writer
  w.sortedAdj = m.adj
  for nbrs in w.sortedAdj.mitems:
    nbrs.sort(proc (x, y: (int, int)): int = cmp(ranks[x[0]], ranks[y[0]]))
  w.inRing = ringBonds(m)
  w.visited = newSeq[bool](n)
  w.bondUsed = newSeq[bool](m.numBonds)
  w.parent = newSeqWith(n, -1)
  w.children = newSeq[seq[(int, int)]](n)
  w.openings = newSeq[seq[int]](n)
  w.closings = newSeq[seq[int]](n)
  w.ringDigit = newSeq[int](m.numBonds)
  let (groups, nGroups) = directionalGroups(m)
  w.dirGroup = groups
  w.groupFlip = newSeqWith(nGroups, -1)
  # components ordered by their lowest-ranked atom
  let comp = connectedComponents(m)
  var starts = newSeqWith(componentCount(comp), -1)
  for i in 0..<n:
    let c = comp[i]
    if starts[c] < 0 or ranks[i] < ranks[starts[c]]: starts[c] = i
  starts.sort(proc (x, y: int): int = cmp(ranks[x], ranks[y]))
  var parts: seq[string]
  for s in starts:
    w.buildTree(s)
    w.output = ""
    w.writeComponent(m, s)
    parts.add w.output
  parts.join(".")

proc toSmiles*(m: Mol; canonical = true; kekule = false): string =
  ## Writes a SMILES string. With `canonical` (default) the result is
  ## independent of the input atom order. With `kekule`, aromatic systems
  ## are written as alternating single/double bonds (uppercase atoms).
  if kekule:
    # kekulize in canonical atom order so the double-bond placement is
    # canonical too
    var k =
      if canonical:
        let ranks = canonicalRanks(m)
        var order = newSeq[int](m.numAtoms)
        for i, r in ranks: order[r] = i
        m.renumberAtoms(order)
      else: m
    kekulize(k)
    toSmilesWithRanks(k, toSeq(0..<k.numAtoms))
  elif canonical:
    toSmilesWithRanks(m, canonicalRanks(m))
  else:
    toSmilesWithRanks(m, toSeq(0..<m.numAtoms))
