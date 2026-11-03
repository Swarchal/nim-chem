import std/[unittest, random, sequtils]
import chem

const corpus = [
  "CCO", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O", "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
  "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "c1ccc2ccccc2c1", "c1ccc2c(c1)ccc1ccccc12",
  "C1CC2CCC1CC2", "C12C3C4C1C5C2C3C45", "OC(=O)C(Cl)(Cl)Cl", "[NH4+].[Cl-]",
  "C[N+](=O)[O-]", "[13CH3]C", "[2H]O[2H]", "c1cc[nH]c1", "O=c1cccc[nH]1",
  "c1ccoc1", "c1ccsc1", "C#N", "N#CC#N", "C=CC=C", "CC(C)(C)C",
  "OC1C(O)C(O)C(O)C(O)C1O", "C1CCC2(CC1)CCCC2", "c1ccc(cc1)-c1ccccc1",
  "CC(C)(C)c1ccc(O)cc1", "CN1CCC[C@H]1c1cccnc1", "C[C@@H](N)C(=O)O",
  "N[C@@H](Cc1ccccc1)C(=O)O", "F/C=C/F", "F/C=C\\F", "CCN(CC)CC",
  "O=S(=O)(O)O", "OP(=O)(O)O", "[Na+].[O-]C(=O)C", "C1=CC=CC=C1",
  "CC1=C(C(=O)CC1)C", "c1ccc2[nH]ccc2c1", "C%10CCCCC%10", "*CC*",
  "CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O", "CC(=O)NC1=CC=C(O)C=C1",
  "[CH3:1][OH:2]", "O.O.O", "C1CC1C1CC1", "c1ccccc1.c1ccccc1",
  "C[S@](=O)c1ccccc1", "CC[C@@](C)(O)CCC", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
  "[Fe+2].[Cl-].[Cl-]", "C1CCCCCCCCCC1",
]

proc shuffled(m: Mol; r: var Rand): Mol =
  var order = toSeq(0..<m.numAtoms)
  r.shuffle(order)
  m.renumberAtoms(order)

suite "canonical SMILES":
  test "alternative spellings agree":
    check canonicalSmiles("OCC") == canonicalSmiles("CCO")
    check canonicalSmiles("C(O)C") == canonicalSmiles("CCO")
    check canonicalSmiles("c1ccccc1O") == canonicalSmiles("Oc1ccccc1")
    check canonicalSmiles("C1CCCCC1C") == canonicalSmiles("CC1CCCCC1")
    check canonicalSmiles("N[C@@H](C)O") == canonicalSmiles("C[C@H](N)O")
    check canonicalSmiles("F\\C=C\\F") == canonicalSmiles("F/C=C/F")
    check canonicalSmiles("Cl.CCO") == canonicalSmiles("OCC.Cl")

  test "stereoisomers differ":
    check canonicalSmiles("C[C@@H](N)O") != canonicalSmiles("C[C@H](N)O")
    check canonicalSmiles("F/C=C\\F") != canonicalSmiles("F/C=C/F")

  test "canonical SMILES is idempotent":
    for s in corpus:
      let c = canonicalSmiles(s)
      checkpoint s & " -> " & c
      check canonicalSmiles(c) == c

  test "invariant under atom renumbering":
    var r = initRand(42)
    for s in corpus:
      let m = parseSmiles(s)
      let c = m.toSmiles
      for _ in 0..<20:
        let p = m.shuffled(r)
        checkpoint s & " canonical " & c
        check p.toSmiles == c
        check p.toSmiles(kekule = true) == m.toSmiles(kekule = true)
        # non-canonical output from a shuffled order must parse back to the same molecule
        check canonicalSmiles(p.toSmiles(canonical = false)) == c

  test "round trip preserves formula":
    for s in corpus:
      let m = parseSmiles(s)
      check parseSmiles(m.toSmiles).molecularFormula == m.molecularFormula
      check parseSmiles(m.toSmiles(canonical = false)).molecularFormula == m.molecularFormula
