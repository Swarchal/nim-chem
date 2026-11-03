import std/[unittest, strutils]
import chem

suite "SMILES parser":
  test "simple chain":
    let m = parseSmiles("CCO")
    check m.numAtoms == 3
    check m.numBonds == 2
    check m.atoms[2].element == 8
    check totalHCount(m, 0) == 3
    check totalHCount(m, 2) == 1

  test "bond orders":
    let m = parseSmiles("C=CC#N")
    check m.bonds[0].order == boDouble
    check m.bonds[1].order == boSingle
    check m.bonds[2].order == boTriple

  test "aromatic ring":
    let m = parseSmiles("c1ccccc1")
    check m.numAtoms == 6
    check m.numBonds == 6
    for bd in m.bonds: check bd.order == boAromatic
    for i in 0..5:
      check m.atoms[i].aromatic
      check totalHCount(m, i) == 1

  test "biaryl bond is single":
    let m = parseSmiles("c1ccccc1c1ccccc1")
    check m.bonds[m.bondBetween(5, 6)].order == boSingle

  test "bracket atoms":
    let m = parseSmiles("[13CH4]")
    check m.atoms[0].isotope == 13
    check m.atoms[0].explicitH == 4
    let n = parseSmiles("[NH4+]")
    check n.atoms[0].charge == 1
    check totalHCount(n, 0) == 4
    let o = parseSmiles("[O--]")
    check o.atoms[0].charge == -2
    let f = parseSmiles("[Fe+3]")
    check f.atoms[0].element == 26
    check f.atoms[0].charge == 3
    let mp = parseSmiles("[CH3:7]C")
    check mp.atoms[0].atomMap == 7
    let se = parseSmiles("c1cc[se]c1")
    check se.atoms[3].element == 34
    check se.atoms[3].aromatic

  test "two letter organic symbols":
    let m = parseSmiles("ClCBr")
    check m.atoms[0].element == 17
    check m.atoms[2].element == 35

  test "branches":
    let m = parseSmiles("CC(C)(C)C")
    check m.degree(1) == 4

  test "ring closures incl %nn":
    let m = parseSmiles("C%12CC%12")
    check m.numBonds == 3
    let r = parseSmiles("C1CC=1")
    check r.bonds[r.bondBetween(0, 2)].order == boDouble

  test "disconnected components":
    let m = parseSmiles("c1ccccc1.Cl")
    check m.numAtoms == 7
    check numComponents(m) == 2

  test "chirality":
    let a = parseSmiles("C[C@H](N)O")
    let b = parseSmiles("C[C@@H](N)O")
    check a.atoms[1].chiral != ctNone
    check a.atoms[1].chiral != b.atoms[1].chiral

  test "wildcard":
    let m = parseSmiles("*C")
    check m.atoms[0].element == 0
    check totalHCount(m, 0) == 0

  test "empty string":
    check parseSmiles("").numAtoms == 0

  test "syntax errors":
    for bad in ["C1CC", "C(C", "CC)", "C==C", "[C", "Xx", "C%1", "(C)",
                "[Zz]", "C=", "C11", "=C"]:
      expect SmilesError:
        discard parseSmiles(bad)

  test "valence errors":
    expect ValenceError:
      discard parseSmiles("C(C)(C)(C)(C)C")
    expect ValenceError:
      discard parseSmiles("O=O=O")
    # sanitize can be disabled
    check parseSmiles("C(C)(C)(C)(C)C", sanitize = false).numAtoms == 6

suite "SMILES parser regressions":
  test "oversized numbers raise SmilesError, not OverflowDefect":
    for bad in ["[99999999999999999999C]", "[C:99999999999999999999]",
                "[CH99999999999999999999]", "[C+99999999999999999999]"]:
      expect SmilesError:
        discard parseSmiles(bad, sanitize = false)

  test "contradictory ring-closure bond directions are rejected":
    expect SmilesError:
      discard parseSmiles("F/C=C/1.F/1")
    # consistent marks (opposite symbols at the two ends) are accepted
    check canonicalSmiles("F/C=C/1.F\\1") == canonicalSmiles("F/C=C/F")

  test "non-tetrahedral chirality gets a clear error":
    for s in ["[C@TB1](F)(Cl)(Br)I", "[C@SP1](F)(Cl)(Br)I", "[C@OH1](F)(Cl)(Br)I"]:
      try:
        discard parseSmiles(s)
        check false
      except SmilesError as e:
        check "only tetrahedral" in e.msg
