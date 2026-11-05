## Molecular descriptors: formula, masses and counts, and the common
## physicochemical descriptors computed as RDKit computes them (Crippen
## logP and molar refractivity, TPSA, Labute's approximate surface area,
## H-bond donors and acceptors, rotatable bonds, ring counts, fraction
## sp3 carbon). See `qed` for drug-likeness.
##
## Each function takes a sanitized molecule (as `parseSmiles` gives).
## Most give the same result whether hydrogens are atoms or counts; like
## RDKit's, rotatable bonds, `numHAcceptors` and Labute ASA change when
## hydrogens are atoms (an O-H is then no longer terminal, and Labute's
## hydrogen term counts H atoms as atoms).

import std/[algorithm, sequtils, tables, math]
import ./[mol, elements, valence, rings, query, smarts, substruct, hydrogens, rdkitdata]

func hydrogenCount*(m: Mol): int =
  ## Total hydrogens: H atoms in the graph plus implicit/explicit H counts.
  for i, a in m.atoms:
    if a.element == 1:
      inc result
    result += totalHCount(m, i)

func numHeavyAtoms*(m: Mol): int =
  for a in m.atoms:
    if a.element != 1:
      inc result

func numAtoms*(m: Mol, includeH: bool): int =
  ## Atom count; with `includeH` implicit and explicit-count H are included.
  result = m.numAtoms
  if includeH:
    for i in 0 ..< m.numAtoms:
      result += totalHCount(m, i)

func totalCharge*(m: Mol): int =
  for a in m.atoms:
    result += a.charge

func molecularFormula*(m: Mol): string =
  ## Molecular formula in Hill order (C, H, then alphabetical; purely
  ## alphabetical if there is no carbon), with a charge suffix such as
  ## "+" or "-2". Isotope labels and wildcard atoms are ignored.
  var counts = initCountTable[string]()
  for a in m.atoms:
    if a.element notin [0, 1]:
      counts.inc symbol(a.element)
  let h = hydrogenCount(m)
  if h > 0:
    counts.inc("H", h)
  func part(sym: string): string =
    let n = counts[sym]
    if n == 1:
      sym
    else:
      sym & $n
  let syms = toSeq(counts.keys).sorted
  if "C" in counts:
    result = part("C")
    if "H" in counts:
      result &= part("H")
    for s in syms:
      if s notin ["C", "H"]:
        result &= part(s)
  else:
    for s in syms:
      result &= part(s)
  result &= chargeSuffix(totalCharge(m))

func massSum(
    m: Mol,
    elementMass: proc(z: int): float {.noSideEffect, nimcall.},
    weighDummies = true,
): float =
  ## Sum of atom masses (labelled isotopes use their exact mass) plus
  ## attached H; without `weighDummies` dummy atoms weigh 0 even with a label.
  for i, a in m.atoms:
    if a.element != 0 or weighDummies:
      result += (
        if a.isotope > 0: isotopeMass(a.element, a.isotope)
        else: elementMass(a.element)
      )
    result += float(totalHCount(m, i)) * elementMass(1)

func rdkitWeight(z: int): float =
  rdkitAtomicWeights[z]

func molWeight*(m: Mol, rdkitWeights = false): float =
  ## Average molecular weight (g/mol), including all hydrogens. With
  ## `rdkitWeights`, RDKit's `MolWt`: its (older) atomic weights, and
  ## dummy atoms weigh nothing.
  if rdkitWeights:
    massSum(m, rdkitWeight, weighDummies = false)
  else:
    massSum(m, avgMass)

func exactMass*(m: Mol): float =
  ## Monoisotopic mass, using the most abundant isotope of each element
  ## (or the labelled isotope). Electron mass is not corrected for charge.
  massSum(m, monoMass)

proc countMatches(m: Mol | PreparedMol, q: Query): int =
  ## Distinct matches (by atom set), unlimited.
  substructMatches(m, q, maxMatches = high(int)).len

# ---------------------------------------------------------------------------
# Wildman-Crippen logP and molar refractivity

var crippen {.threadvar.}:
  tuple[queries: seq[Query], elements: seq[set[uint8]], needs: set[QueryNeeds]]

proc crippenTyping(m: Mol): (Mol, seq[int]) =
  ## `m` with all H as atoms, and each atom's index in `crippenTypes` (-1
  ## if none fits): the first type whose pattern matches with its first
  ## atom there.
  if crippen.queries.len == 0:
    crippen.queries.parsedOnce(crippenTypes.mapIt(it.smarts))
    for q in crippen.queries:
      crippen.elements.add elements(q.atoms[0].expr)
      crippen.needs = crippen.needs + q.needs
  let p = prepare(addHs(m), crippen.needs)
  let n = p.mol.numAtoms
  var types = newSeqWith(n, -1)
  var left = n
  for k, q in crippen.queries:
    if left == 0:
      break
    # only untyped atoms of an element the type allows
    var candidates: seq[int]
    for a in 0 ..< n:
      if types[a] < 0 and uint8(p.mol.atoms[a].element) in crippen.elements[k]:
        candidates.add a
    if candidates.len == 0:
      continue
    for j, hit in p.matchesAt(q, candidates):
      if hit:
        types[candidates[j]] = k
        dec left
  (p.mol, types)

proc crippenAtomTypes*(m: Mol): seq[string] =
  ## Each atom's Wildman-Crippen type label ("C1", "N11", ...; "" if none
  ## fits), as RDKit assigns them. Hydrogens on an atom are typed too (they
  ## add to its contributions) but aren't listed.
  let (_, types) = crippenTyping(m)
  result = newSeq[string](m.numAtoms)
  for i in 0 ..< m.numAtoms:
    if types[i] >= 0:
      result[i] = crippenTypes[types[i]].label

proc crippenContribs*(m: Mol): seq[tuple[logP, mr: float]] =
  ## Per-atom contributions to `crippenLogP` and `crippenMR`, with each
  ## hydrogen's contribution included in its neighbour's (an H atom of `m`
  ## keeps its own).
  let (h, types) = crippenTyping(m)
  result = newSeq[tuple[logP, mr: float]](m.numAtoms)
  for a, k in types:
    if k < 0:
      continue
    # added H atoms come after `m`'s atoms and have one neighbour
    let owner =
      if a < m.numAtoms:
        a
      else:
        h.adj[a][0][0]
    result[owner].logP += crippenTypes[k].logP
    result[owner].mr += crippenTypes[k].mr

proc crippenLogP*(m: Mol): float =
  ## Wildman-Crippen octanol/water partition coefficient (RDKit's
  ## `MolLogP`): a sum of atom-type contributions over all atoms,
  ## hydrogens included. (Sum `crippenContribs` for logP and MR at once.)
  crippenContribs(m).mapIt(it.logP).sum

proc crippenMR*(m: Mol): float =
  ## Wildman-Crippen molar refractivity (RDKit's `MolMR`).
  crippenContribs(m).mapIt(it.mr).sum

# ---------------------------------------------------------------------------
# Polar surface area

func inThreeRing(m: Mol, i: int): bool =
  for x in m.neighbors(i):
    for y in m.neighbors(i):
      if x < y and m.bondBetween(x, y) >= 0:
        return true

func tpsaContribs*(m: Mol): seq[float] =
  ## Per-atom contributions to `tpsa` (nonzero only for N and O).
  result = newSeq[float](m.numAtoms)
  for i, a in m.atoms:
    if a.element notin [7, 8]:
      continue
    let hs = totalHCount(m, i, includeNeighbors = true)
    var nbrs, nSing, nDoub, nTrip, nArom = 0
    for (nb, bi) in m.atomBonds(i):
      if m.atoms[nb].element == 1:
        continue
      inc nbrs
      case m.bonds[bi].order
      of boSingle:
        inc nSing
      of boDouble:
        inc nDoub
      of boTriple:
        inc nTrip
      of boAromatic:
        inc nArom
      else:
        discard
    let q = a.charge
    template rule(cond: bool, value: float) =
      if cond:
        return value

    proc nitrogen(): float =
      case nbrs
      of 1:
        rule(hs == 0 and nTrip == 1 and q == 0, 23.79)
        rule(hs == 1 and nDoub == 1 and q == 0, 23.85)
        rule(hs == 2 and nSing == 1 and q == 0, 26.02)
        rule(hs == 2 and nDoub == 1 and q == 1, 25.59)
        rule(hs == 3 and nSing == 1 and q == 1, 27.64)
      of 2:
        rule(hs == 0 and nSing == 1 and nDoub == 1 and q == 0, 12.36)
        rule(hs == 0 and nTrip == 1 and nDoub == 1 and q == 0, 13.60)
        rule(
          hs == 1 and nSing == 2 and q == 0, (if inThreeRing(m, i): 21.94 else: 12.03)
        )
        rule(hs == 0 and nTrip == 1 and nSing == 1 and q == 1, 4.36)
        rule(hs == 1 and nDoub == 1 and nSing == 1 and q == 1, 13.97)
        rule(hs == 2 and nSing == 2 and q == 1, 16.61)
        rule(hs == 0 and nArom == 2 and q == 0, 12.89)
        rule(hs == 1 and nArom == 2 and q == 0, 15.79)
        rule(hs == 1 and nArom == 2 and q == 1, 14.14)
      of 3:
        rule(hs == 0 and nSing == 3 and q == 0, (if inThreeRing(m, i): 3.01 else: 3.24))
        rule(hs == 0 and nSing == 1 and nDoub == 2 and q == 0, 11.68)
        rule(hs == 0 and nSing == 2 and nDoub == 1 and q == 1, 3.01)
        rule(hs == 1 and nSing == 3 and q == 1, 4.44)
        rule(hs == 0 and nArom == 3 and q == 0, 4.41)
        rule(hs == 0 and nSing == 1 and nArom == 2 and q == 0, 4.93)
        rule(hs == 0 and nDoub == 1 and nArom == 2 and q == 0, 8.39)
        rule(hs == 0 and nArom == 3 and q == 1, 4.10)
        rule(hs == 0 and nSing == 1 and nArom == 2 and q == 1, 3.88)
      of 4:
        rule(hs == 0 and nSing == 4 and q == 1, 0.0)
      else:
        discard
      max(30.5 - float(nbrs) * 8.2 + float(hs) * 1.5, 0.0)

    proc oxygen(): float =
      case nbrs
      of 1:
        rule(hs == 0 and nDoub == 1 and q == 0, 17.07)
        rule(hs == 1 and nSing == 1 and q == 0, 20.23)
        rule(hs == 0 and nSing == 1 and q == -1, 23.06)
      of 2:
        rule(
          hs == 0 and nSing == 2 and q == 0, (if inThreeRing(m, i): 12.53 else: 9.23)
        )
        rule(hs == 0 and nArom == 2 and q == 0, 13.14)
      else:
        discard
      max(28.5 - float(nbrs) * 8.6 + float(hs) * 1.5, 0.0)

    result[i] =
      if a.element == 7:
        nitrogen()
      else:
        oxygen()

func tpsa*(m: Mol): float =
  ## Topological polar surface area (Ertl et al. 2000), from N and O only,
  ## as RDKit's `TPSA` and Daylight's tpsa.c compute it.
  sum(tpsaContribs(m))

func labuteASA*(m: Mol): float =
  ## Labute's approximate surface area (RDKit's `LabuteASA`). As RDKit
  ## computes it, each atom is reduced by one hydrogen's overlap and one
  ## hydrogen sphere is added, whatever the hydrogen count (none if all
  ## atoms are dummies, whose radius is 0).
  let rH = rb0Radii[1]
  let rads = m.atoms.mapIt(rb0Radii[it.element])
  var v = newSeq[float](m.numAtoms)
  func overlap(ri, rj, bij: float): (float, float) =
    ## how much of each sphere a bond of length `bij` buries
    let dij = min(max(abs(ri - rj), bij), ri + rj)
    if dij == 0:
      return (0.0, 0.0) # two dummy atoms
    (rj * rj - (ri - dij) ^ 2 / dij, ri * ri - (rj - dij) ^ 2 / dij)
  for bd in m.bonds:
    let (ri, rj) = (rads[bd.a], rads[bd.b])
    let shorter =
      case bd.order
      of boAromatic: 0.1
      of boDouble: 0.2
      of boTriple: 0.3
      else: 0.0
    let (vi, vj) = overlap(ri, rj, ri + rj - shorter)
    v[bd.a] += vi
    v[bd.b] += vj
  var vh = 0.0
  for i, ri in rads:
    let (vi, vj) = overlap(ri, rH, ri + rH)
    v[i] += vi
    vh += vj
  for i, ri in rads:
    result += 4 * PI * ri * ri - PI * ri * v[i]
  if vh != 0:
    result += 4 * PI * rH * rH - PI * rH * vh

# ---------------------------------------------------------------------------
# Lipinski-style counts

type CountPattern = enum
  cpDonor
  cpAcceptor
  cpRotatable
  cpStrictRotatable

const
  donorSmarts = "[$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])]"
  acceptorSmarts =
    "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-])," &
    "$([N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0,o,s;+0;X2])]"
  rotatable = "[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]"
  # RDKit's strict definition: also not bonds to CX3 (CF3, CCl3, CBr3, tBu)
  # or to the C of an amide-like group, nor a methyl's bond
  strictEnd =
    "!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&" &
    "!$(C([CH3])([CH3])[CH3])&!$([CH3])"
  strictRotatable =
    "[" & strictEnd & "&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&" &
    "!$([#7,O,S!D1]-!@[CD3]=[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&" &
    "!$([#7!D1]-!@[CD3]=[N+])]-,:;!@[" & strictEnd & "]"

  countSmarts: array[CountPattern, string] =
    [donorSmarts, acceptorSmarts, rotatable, strictRotatable]

var countQueries {.threadvar.}: seq[Query]

proc countQuery(k: CountPattern): lent Query =
  countQueries.parsedOnce(countSmarts)
  countQueries[ord(k)]

proc numHDonors*(m: Mol | PreparedMol): int =
  ## H-bond donors, by RDKit's definition (`NumHDonors`): N with H (and
  ## valence 3, or 4 if +1), neutral O and S with one H, aromatic nH.
  countMatches(m, countQuery(cpDonor))

proc numHAcceptors*(m: Mol | PreparedMol): int =
  ## H-bond acceptors, by RDKit's definition (`NumHAcceptors`): O and S
  ## (hydroxyl, thiol, ether, anions; not OH on a C=O), N other than amide
  ## N, two-connected aromatic n, aromatic o and s.
  countMatches(m, countQuery(cpAcceptor))

func numLipinskiHDonors*(m: Mol): int =
  ## Hydrogens on N and O (`NHOHCount`, the rule-of-five donor count).
  for i, a in m.atoms:
    if a.element in [7, 8]:
      result += totalHCount(m, i, includeNeighbors = true)

func numLipinskiHAcceptors*(m: Mol): int =
  ## N and O atoms (`NOCount`, the rule-of-five acceptor count).
  m.atoms.countIt(it.element in [7, 8])

func numHeteroatoms*(m: Mol): int =
  ## Atoms other than C and H.
  m.atoms.countIt(it.element notin [1, 6])

proc numRotatableBonds*(m: Mol | PreparedMol, strict = true): int =
  ## Rotatable bonds, as RDKit counts them: single, non-ring bonds
  ## between two non-terminal atoms, neither in a triple bond. `strict`
  ## (RDKit's default) also excludes bonds to CF3, CCl3, CBr3 and t-butyl
  ## groups and the C-N bond of amides and similar.
  countMatches(m, countQuery(if strict: cpStrictRotatable else: cpRotatable))

proc lipinskiViolations*(m: Mol): int =
  ## How many of Lipinski's rule-of-five limits `m` exceeds: molecular
  ## weight 500, logP 5 (Crippen), 5 donors and 10 acceptors (the NHOH
  ## and NO counts).
  int(m.molWeight > 500) + int(m.crippenLogP > 5) + int(m.numLipinskiHDonors > 5) +
    int(m.numLipinskiHAcceptors > 10)

# ---------------------------------------------------------------------------
# Rings and carbon saturation

type RingCounts* = object
  ## Counts over the symmetrized SSSR (see `rings.symmetrizedSssr`), as
  ## RDKit's ring descriptors count them. A ring is aromatic if all its
  ## bonds are, saturated if all are single, a heterocycle if any atom
  ## isn't carbon; aliphatic means not aromatic.
  rings*, aromatic*, aliphatic*, saturated*, heterocycles*: int
  aromaticCarbocycles*, aromaticHeterocycles*: int
  aliphaticCarbocycles*, aliphaticHeterocycles*: int
  saturatedCarbocycles*, saturatedHeterocycles*: int

proc ringCounts*(m: Mol): RingCounts =
  for ring in symmetrizedSssr(m):
    let bonds = m.ringBondIndices(ring)
    let arom = bonds.allIt(m.bonds[it].order == boAromatic)
    let sat = bonds.allIt(m.bonds[it].order == boSingle)
    let hetero = ring.anyIt(m.atoms[it].element != 6)
    template tally(total, heteroCount, carboCount: var int) =
      inc total
      if hetero:
        inc heteroCount
      else:
        inc carboCount

    inc result.rings
    if hetero:
      inc result.heterocycles
    if arom:
      tally(result.aromatic, result.aromaticHeterocycles, result.aromaticCarbocycles)
    else:
      tally(result.aliphatic, result.aliphaticHeterocycles, result.aliphaticCarbocycles)
    if sat:
      tally(result.saturated, result.saturatedHeterocycles, result.saturatedCarbocycles)

func fractionCSP3*(m: Mol): float =
  ## Fraction of carbon atoms that are sp3 (0 without carbon), taken as
  ## RDKit does: carbon with four neighbours, hydrogens included. (So a
  ## carbanion R3C- doesn't count, nor a radical.)
  var c, sp3 = 0
  for i, a in m.atoms:
    if a.element == 6:
      inc c
      if m.degree(i) + totalHCount(m, i) == 4:
        inc sp3
  if c == 0:
    0.0
  else:
    sp3 / c
