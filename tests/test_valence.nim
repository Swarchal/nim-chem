import std/unittest
import chem

proc hs(s: string): seq[int] =
  let m = parseSmiles(s)
  for i in 0..<m.numAtoms: result.add totalHCount(m, i)

suite "valence":
  test "organic subset":
    check hs("C") == @[4]
    check hs("N") == @[3]
    check hs("O") == @[2]
    check hs("Cl") == @[1]
    check hs("C=O") == @[2, 0]
    check hs("C#N") == @[1, 0]

  test "hypervalent defaults":
    check hs("CS(=O)(=O)C") == @[3, 0, 0, 0, 3]
    check hs("OP(=O)(O)O") == @[1, 0, 0, 1, 1]
    check hs("CN(=O)=O") == @[3, 0, 0, 0]

  test "aromatic":
    check hs("c1ccccc1") == @[1, 1, 1, 1, 1, 1]
    check hs("c1ccncc1") == @[1, 1, 1, 0, 1, 1]
    check hs("c1ccoc1") == @[1, 1, 1, 0, 1]
    check hs("c1ccsc1") == @[1, 1, 1, 0, 1]
    check hs("c1cc[nH]c1") == @[1, 1, 1, 1, 1]
    check hs("O=c1cccc[nH]1") == @[0, 0, 1, 1, 1, 1, 1]

  test "bracket atoms have no implicit H":
    check hs("[C]") == @[0]
    check hs("[CH3]") == @[3]

  test "charged atoms":
    check validateValences(parseSmiles("C[N+](C)(C)C")).len == 0
    check validateValences(parseSmiles("C[N+](=O)[O-]")).len == 0
    check validateValences(parseSmiles("[N+](C)(C)(C)(C)C", sanitize = false)).len == 1

  test "metals are not validated":
    check validateValences(parseSmiles("[Fe](C)(C)(C)(C)(C)C")).len == 0
