## Morgan (circular) fingerprints, identical to RDKit's: the same bit IDs
## as `GetMorganGenerator` / the legacy `GetMorganFingerprint`, so
## fingerprints can be compared with ones RDKit computed.
##
## Each atom starts with an invariant (from its element, degree, H count,
## charge, isotope and ring membership; or, for feature-based "FCFP"
## fingerprints, its pharmacophore roles). Each round hashes in the
## neighbours' invariants, so round `r` describes the atom's environment out
## to `r` bonds. Every distinct environment contributes its hash as a bit
## ID; an environment covering the same bonds as one seen before is
## dropped. Radius 2 corresponds to ECFP4 / FCFP4, radius 3 to ECFP6.
##
## Chirality is not included (RDKit's `useChirality` uses CIP labels).

import std/[tables, sets, algorithm, bitops]
import ./[mol, elements, valence, rings, query, smarts, substruct, bitvect, rdkitdata]

type MorganEnv* = tuple
  ## One environment: its bit ID (before folding), centre atom and radius.
  id: uint32
  atom, radius: int

func hashCombine(seed: var uint32; v: uint32) {.inline.} =
  ## boost::hash_combine, as RDKit uses it (on 32-bit values).
  seed = seed xor (v + 0x9e3779b9'u32 + (seed shl 6) + (seed shr 2))

func rdkitBondType(o: BondOrder): uint32 =
  ## RDKit's `Bond::BondType` numbering.
  case o
  of boNone: 0
  of boSingle: 1
  of boDouble: 2
  of boTriple: 3
  of boQuadruple: 4
  of boAromatic: 12

func massShift(a: Atom): int =
  ## RDKit's whole-number mass difference of an atom from its element's
  ## average weight: 0 without an isotope label, and for dummy atoms (whose
  ## mass RDKit takes as 0, label or not). The shift is truncated, so a few
  ## isotopes (32P, 28Al, ...) need RDKit's own weights to match it.
  if a.isotope == 0 or a.element == 0: 0
  else: int(isotopeMass(a.element, a.isotope) - rdkitAtomicWeights[a.element])

proc connectivityInvariants*(m: Mol): seq[uint32] =
  ## Per-atom starting invariants for ECFP-style fingerprints (RDKit's
  ## `GetConnectivityInvariants`): atomic number, degree including H, H
  ## count, charge, isotope mass shift and ring membership.
  let inRing = ringAtoms(m)
  result = newSeq[uint32](m.numAtoms)
  for i, a in m.atoms:
    let h = totalHCount(m, i)
    var hAtoms = 0
    for nb in m.neighbors(i):
      if m.atoms[nb].element == 1: inc hAtoms
    var x = 0'u32
    for v in [a.element.uint32, uint32(m.degree(i) + h), uint32(h + hAtoms),
              cast[uint32](a.charge.int32), cast[uint32](a.massShift.int32)]:
      hashCombine(x, v)
    if inRing[i]: hashCombine(x, 1)
    result[i] = x

const featureSmarts = [
  # donor
  "[$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),n&H1&+0]",
  # acceptor
  "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-])," &
    "$([N;v3;!$(N-*=[O,N,P,S])]),n&H0&+0,$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]",
  # aromatic
  "[a]",
  # halogen
  "[F,Cl,Br,I]",
  # basic
  "[#7;+,$([N;H2&+0][$([C,a]);!$([C,a](=O))])," &
    "$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))])," &
    "$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]",
  # acidic
  "[$([C,S](=[O,S,P])-[O;H1,-1])]"]

var featurePatterns {.threadvar.}: seq[Query]   # parsed on first use

proc featureInvariants*(m: Mol): seq[uint32] =
  ## Per-atom starting invariants for FCFP-style fingerprints (RDKit's
  ## `GetFeatureInvariants`): a bit set for each role the atom plays, in
  ## the order donor, acceptor, aromatic, halogen, basic, acidic. As in
  ## RDKit, only the first 1,000 atoms matching a pattern get its bit.
  featurePatterns.parsedOnce(featureSmarts)
  result = newSeq[uint32](m.numAtoms)
  for j, q in featurePatterns:
    # RDKit's SubstructMatch default: at most 1,000 matches per pattern
    for match in m.substructMatches(q, uniquify = false, maxMatches = 1000):
      result[match[0]] = result[match[0]] or (1'u32 shl j)

func morganEnvironments*(m: Mol; invariants: openArray[uint32]; radius = 2;
                         useBondTypes = true): seq[MorganEnv] =
  ## The environments of a Morgan fingerprint built from per-atom
  ## `invariants`: every atom at radius 0, then each round's new
  ## environments (the same set as RDKit's, though within a round RDKit
  ## orders them differently).
  let n = m.numAtoms
  if invariants.len != n:
    raise newException(ValueError, "need one invariant per atom")
  # per atom, the bonds within the current radius: `words` words each
  let words = (m.numBonds + 63) div 64
  var current = @invariants
  for i in 0..<n: result.add (current[i], i, 0)
  var hood, nextHood = newSeq[uint64](n * words)
  var seen = initHashSet[seq[uint64]]()
  var key = newSeq[uint64](words)
  var dead = newSeq[bool](n)
  var next = newSeq[uint32](n)
  var nbrs: seq[(uint32, uint32)]
  var round: seq[(uint32, int)]
  for layer in 0..<radius:
    next.fill(0)   # dead atoms pass on 0, as in RDKit
    nextHood = hood
    round.setLen 0
    for i in 0..<n:
      if dead[i]: continue
      if m.degree(i) == 0:
        dead[i] = true
        continue
      nbrs.setLen 0
      for (nb, bi) in m.adj[i]:
        nextHood[i * words + (bi shr 6)].setBit(bi and 63)
        for k in 0..<words: nextHood[i * words + k] = nextHood[i * words + k] or hood[nb * words + k]
        let bt = if useBondTypes: rdkitBondType(m.bonds[bi].order) else: 1
        nbrs.add (bt, current[nb])
      nbrs.sort()
      var h = layer.uint32
      hashCombine(h, current[i])
      for (bt, inv) in nbrs:
        var pairHash = 0'u32
        hashCombine(pairHash, bt)
        hashCombine(pairHash, inv)
        hashCombine(h, pairHash)
      next[i] = h
      round.add (h, i)
    # an environment covering the same bonds as an earlier one adds nothing
    # (of equal ones this round, the smallest hash is kept), and its atom
    # takes no further part
    round.sort()
    for (h, i) in round:
      for k in 0..<words: key[k] = nextHood[i * words + k]
      if seen.containsOrIncl(key): dead[i] = true
      else: result.add (h, i, layer + 1)
    swap(current, next)
    swap(hood, nextHood)

proc morganInvariants*(m: Mol; useFeatures = false): seq[uint32] =
  ## The starting invariants Morgan fingerprints use: `featureInvariants`
  ## or `connectivityInvariants`.
  if useFeatures: featureInvariants(m) else: connectivityInvariants(m)

proc morganCounts*(m: Mol; radius = 2; useFeatures = false;
                   useBondTypes = true): CountTable[uint32] =
  ## The unfolded Morgan fingerprint: how often each environment hash
  ## occurs (RDKit's `GetSparseCountFingerprint`).
  for e in morganEnvironments(m, morganInvariants(m, useFeatures), radius, useBondTypes):
    result.inc e.id

proc morganFingerprint*(m: Mol; radius = 2; nBits = 2048; useFeatures = false;
                        useBondTypes = true): BitVect =
  ## The Morgan fingerprint folded to `nBits` bits (RDKit's
  ## `GetFingerprint`; radius 2 is ECFP4, or FCFP4 with `useFeatures`).
  if nBits <= 0: raise newException(ValueError, "nBits must be positive")
  result = initBitVect(nBits)
  for e in morganEnvironments(m, morganInvariants(m, useFeatures), radius, useBondTypes):
    result.incl int(e.id mod nBits.uint32)
