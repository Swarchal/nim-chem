import std/[unittest, strutils]
import chem

suite "SMILES writer":
  test "non-canonical follows input order":
    check parseSmiles("OCC").toSmiles(canonical = false) == "OCC"
    check parseSmiles("C(O)C").toSmiles(canonical = false) == "C(O)C"

  test "atom tokens":
    check canonicalSmiles("[NH4+]") == "[NH4+]"
    check canonicalSmiles("[O--]") == "[O-2]"
    check canonicalSmiles("[13CH4]") == "[13CH4]"
    check canonicalSmiles("[CH3]") == "[CH3]"
    check canonicalSmiles("[CH4]") == "C"
    check canonicalSmiles("[nH]1cccc1") == "c1cc[nH]c1"
    check canonicalSmiles("[CH3:1]O") == "[CH3:1]O"

  test "bond tokens":
    check canonicalSmiles("C=C") == "C=C"
    check canonicalSmiles("C#N") == "C#N"
    check canonicalSmiles("c1ccccc1-c1ccccc1").count('-') == 1

  test "many ring closures use %nn":
    # a ladder-like molecule with 12 simultaneously open rings
    let m = parseSmiles("C12C3C4C5C6C7C8C9C%10C%11C%12CC%12C%11C%10C9C8C7C6C5C4C3C1C2", sanitize = false)
    let out1 = m.toSmiles(canonical = false)
    check "%1" in out1
    check parseSmiles(out1, sanitize = false).molecularFormula == m.molecularFormula
    check canonicalSmiles(out1) == canonicalSmiles(m.toSmiles)

  test "disconnected":
    check canonicalSmiles("[Na+].[Cl-]") == canonicalSmiles("[Cl-].[Na+]")

suite "SMILES writer regressions":
  test "non-ring aromatic bond round-trips":
    let m = parseSmiles("c1ccccc1:c1ccccc1")
    let out1 = m.toSmiles
    check ":" in out1
    let back = parseSmiles(out1)
    var nArom = 0
    for bd in back.bonds:
      if bd.order == boAromatic: inc nArom
    check nArom == 13
    check out1 != canonicalSmiles("c1ccccc1-c1ccccc1")

  test "unspecified bond order raises instead of writing '~'":
    var m: Mol
    m.addAtom(initAtom(6))
    m.addAtom(initAtom(6))
    m.addBond(0, 1, boNone)
    expect ChemError:
      discard m.toSmiles

  test "very long chain does not overflow the stack":
    const n = 300_000
    let m = parseSmiles("C".repeat(n))
    let s = m.toSmiles(canonical = false)
    check s.len == n
    check parseSmiles(s).numAtoms == n
