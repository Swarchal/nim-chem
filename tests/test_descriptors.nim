import std/[unittest, math]
import chem

suite "descriptors":
  test "formula":
    check parseSmiles("CC(=O)Oc1ccccc1C(=O)O").molecularFormula == "C9H8O4"
    check parseSmiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C").molecularFormula == "C8H10N4O2"
    check parseSmiles("[NH4+]").molecularFormula == "H4N+"
    check parseSmiles("O").molecularFormula == "H2O"
    check parseSmiles("[O-]S(=O)(=O)[O-]").molecularFormula == "O4S-2"
    check parseSmiles("[2H]C([2H])([2H])Cl").molecularFormula == "CH3Cl"

  test "molecular weight":
    check abs(parseSmiles("CC(=O)Oc1ccccc1C(=O)O").molWeight - 180.159) < 1e-3
    check abs(parseSmiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C").molWeight - 194.194) < 1e-3
    check abs(parseSmiles("O").molWeight - 18.015) < 1e-3

  test "exact mass":
    check abs(parseSmiles("O").exactMass - 18.010565) < 1e-5
    check abs(parseSmiles("CC(=O)Oc1ccccc1C(=O)O").exactMass - 180.042259) < 1e-5
    check abs(parseSmiles("[13CH4]").exactMass - 17.034655) < 1e-5

  test "counts":
    let m = parseSmiles("CC(=O)Oc1ccccc1C(=O)O")
    check m.numHeavyAtoms == 13
    check m.numAtoms(includeH = true) == 21
    check m.numBonds == 13
    check m.numRings == 1
    check parseSmiles("c1ccc2ccccc2c1").numRings == 2
    check parseSmiles("C1CC2CCC1CC2").numRings == 2
    check parseSmiles("CCO.CC").numRings == 0

  test "rotatable bonds":
    check parseSmiles("CCCC").numRotatableBonds == 1
    check parseSmiles("c1ccccc1-c1ccccc1").numRotatableBonds == 1
    check parseSmiles("C1CCCCC1").numRotatableBonds == 0
