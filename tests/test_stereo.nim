import std/[unittest, random, sequtils, sets]
import chem
import ./helpers

proc variants(s: string; n = 50): HashSet[string] =
  ## Canonical SMILES of `s` under random atom orders, written both directly
  ## and via a non-canonical round trip.
  var r = initRand(11)
  let m = parseSmiles(s)
  for _ in 0..<n:
    let p = m.shuffled(r)
    result.incl p.toSmiles
    result.incl canonicalSmiles(p.toSmiles(canonical = false))

proc numChiral(s: string): int =
  for a in parseSmiles(s).atoms:
    if a.chiral != ctNone: inc result

proc numStereoBonds(s: string): int =
  for bd in parseSmiles(s).bonds:
    if bd.stereo != bsNone: inc result

const unique = [
  "C[C@H](O)[C@@H](C)O",                       # meso-2,3-butanediol
  "C[C@H](O)[C@H](C)O", "C[C@@H](O)[C@@H](C)O",
  "OC(=O)[C@H](O)[C@@H](O)C(=O)O",             # (R,R)-tartaric acid
  "OC(=O)[C@@H](O)[C@H](O)C(=O)O",             # (S,S)
  "OC(=O)[C@H](O)[C@H](O)C(=O)O",              # meso
  "C[C@@H]1CC[C@H](C)CC1", "C[C@@H]1CC[C@@H](C)CC1",   # cis/trans ring
  "OC(=O)[C@H](O)[C@H](O)[C@H](O)C(=O)O",      # pseudoasymmetric centre
  "C/C=C/C=C/C", "C/C=C\\C=C/C", "F/C(Cl)=C/F", "F/C=C\\F",
  "C[S@](=O)c1ccccc1", "N[C@@H](Cc1ccccc1)C(=O)O",
  "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
  "C1/C=C\\CCCCC1", "C/C=C/c1ccccc1",
  # linked rings whose ring stereo is only labelled after a later tie-break
  "C[C@H]1CC[C@H]([C@H]2CC[C@@H](C)CC2)CC1",
  "C[C@H]1CC[C@H](C[C@H]2CC[C@@H](C)CC2)CC1",
  "C[C@H]1C[C@H]([C@H]2C[C@@H](C)C2)C1",
]

suite "stereo canonical SMILES":
  test "unique under atom renumbering":
    for s in unique:
      checkpoint s
      check variants(s).len == 1

  test "stereoisomers stay distinct":
    # (R,R), (S,S) and meso tartaric acid (CIP labels checked with RDKit)
    let tartaric = ["OC(=O)[C@H](O)[C@@H](O)C(=O)O",
                    "OC(=O)[C@@H](O)[C@H](O)C(=O)O",
                    "OC(=O)[C@H](O)[C@H](O)C(=O)O"].mapIt(canonicalSmiles(it))
    check tartaric.toHashSet.len == 3
    # both spellings of meso are the same molecule
    check canonicalSmiles("OC(=O)[C@H](O)[C@H](O)C(=O)O") ==
          canonicalSmiles("OC(=O)[C@@H](O)[C@@H](O)C(=O)O")
    check canonicalSmiles("C[C@@H]1CC[C@H](C)CC1") != canonicalSmiles("C[C@@H]1CC[C@@H](C)CC1")
    check canonicalSmiles("F/C=C/F") != canonicalSmiles("F/C=C\\F")
    check canonicalSmiles("C/C=C/C=C/C") != canonicalSmiles("C/C=C\\C=C/C")

  test "equivalent double-bond spellings agree":
    let trans = canonicalSmiles("F/C=C/F")
    for s in ["F\\C=C\\F", "C(\\F)=C/F", "F/C=C/1.F\\1"]:
      checkpoint s
      check canonicalSmiles(s) == trans
    let a = canonicalSmiles("F/C(Cl)=C/F")
    for s in ["Cl\\C(F)=C/F", "F/C(/Cl)=C/F", "F/C=C(/F)Cl"]:
      checkpoint s
      check canonicalSmiles(s) == a

  test "equivalent tetrahedral spellings agree":
    check canonicalSmiles("N[C@@H](C)C(=O)O") == canonicalSmiles("C[C@H](N)C(=O)O")
    check canonicalSmiles("C[S@](=O)c1ccccc1") == canonicalSmiles("[S@](C)(=O)c1ccccc1")
    check canonicalSmiles("C[S@](=O)c1ccccc1") != canonicalSmiles("C[S@@](=O)c1ccccc1")

suite "stereo perception":
  test "meaningless tags are removed":
    check canonicalSmiles("CC[C@@](C)(O)CC") == "CCC(C)(O)CC"
    check canonicalSmiles("F/C=C(/F)F") == "FC=C(F)F"
    check numChiral("C[N@](CC)CCC") == 0         # amine nitrogen inverts
    check numChiral("C[C@](=O)O") == 0           # sp2 carbon
    check numStereoBonds("C/C=C") == 0           # marks on one end only
    check numStereoBonds("C/C=C(C)/C") == 0      # identical substituents
    check numStereoBonds("C1/C=C\\CCCC1") == 0   # double bond in a 6-ring

  test "real stereo is kept":
    check numChiral("C[S@](=O)c1ccccc1") == 1
    check numChiral("C[C@@H]1CC[C@H](C)CC1") == 2
    check numStereoBonds("C1/C=C\\CCCCC1") == 1  # 8-ring
    check numStereoBonds("C/C=C/C=C/C") == 2

  test "pseudoasymmetric centre":
    # 2,3,4-trihydroxyglutaric acid: the middle carbon is stereogenic only
    # when its two neighbours have opposite configurations
    check numChiral("OC(=O)[C@H](O)[C@H](O)[C@H](O)C(=O)O") == 3
    check numChiral("OC(=O)[C@H](O)[C@H](O)[C@@H](O)C(=O)O") == 2

  test "an explicit [H] counts as the implicit H":
    check numChiral("[H][C@H](F)Cl") == 0
    check numChiral("C[C@@H]([H])O") == 0
    check canonicalSmiles("[H][C@H](F)Cl") == canonicalSmiles("[H][C@@H](F)Cl")
    check numChiral("[H][C@](F)(Cl)Br") == 1
    check numChiral("[2H][C@H](F)Cl") == 1      # deuterium is distinct

  test "contradictory marks raise SmilesError":
    expect SmilesError:
      discard parseSmiles("C(/F)(/Cl)=C/F")

  test "stereo errors point at the offending input":
    try:
      discard parseSmiles("CCCCCC(/F)(/Cl)=C/F")
      fail()
    except SmilesError as e:
      check e.pos in 7..11          # one of the two conflicting marks

  test "subMol rejects bad atom lists":
    let m = parseSmiles("CC(F)Cl")
    expect ChemError: discard m.subMol([0, 1, 1])
    expect ChemError: discard m.subMol([0, 7])
    check m.subMol([1, 2]).toSmiles == "CF"

  test "renumbering keeps double-bond stereo":
    let m = parseSmiles("F/C=C/Cl")
    let p = m.renumberAtoms([3, 2, 1, 0])
    check p.toSmiles == m.toSmiles

suite "three-coordinate centres and bridgeheads":
  test "aziridine N spellings agree (RDKit ring-digit convention)":
    # all of these are RDKit's own random SMILES for C[N@]1CC1(C)C
    let want = canonicalSmiles("C[N@]1CC1(C)C")
    for s in ["C1[N@](C1(C)C)C", "CC1(C)C[N@@]1C", "[N@]1(C)C(C)(C)C1",
              "C1C(C)(C)[N@]1C", "[N@@]1(C)CC1(C)C", "C1[N@@](C)C1(C)C"]:
      checkpoint s
      check canonicalSmiles(s) == want
    check canonicalSmiles("C[N@@]1CC1(C)C") != want

  test "sulfoxides follow OpenSMILES, also in rings":
    # expected value is RDKit's canonical SMILES for the same input
    check canonicalSmiles("O=[S@]1CCC[C@H]1C") == "C[C@@H]1CCC[S@@]1=O"
    check variants("O=[S@]1CCC[C@H]1C", 20).len == 1

  test "bridgehead stereo":
    # brucine: the bridgehead N is kept, and its tag is read consistently
    let brucine = "COc1cc2c(cc1OC)[C@]13CC[N@@]4CC5=CCO[C@H]6CC(=O)N2[C@@H]1[C@H]6[C@@H]5[C@H]3C4"
    check variants(brucine, 20).len == 1
    var nChiralN = 0
    for a in parseSmiles(brucine).atoms:
      if a.chiral != ctNone and a.element == 7: inc nChiralN
    check nChiralN == 1
    # bicyclo[2.2.2]octane bridgehead carbons are kept (ring stereo)
    check numChiral("CCOC(=O)[C@@H]1[C@H]2CC[C@@H](CC2)[C@H]1C(=O)OCC") == 4
    # cage aminals: every N is fixed by the cage, so the tags are dropped
    check numChiral("C1C[N@]2C[N@@]3CC[N@](C[N@]1C2)C3") == 0
