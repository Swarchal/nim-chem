import std/[unittest, tables, sequtils, algorithm, math, random, strutils]
import chem
import ./helpers

proc ids(c: CountTable[uint32]): seq[(uint32, int)] =
  result = toSeq(c.pairs)
  result.sort()

suite "bit vectors":
  test "setting, counting and iterating":
    var v = initBitVect(130)
    check v.len == 130 and v.count == 0
    for i in [0, 63, 64, 129]: v.incl i
    check v.count == 4
    check toSeq(v.onBits) == @[0, 63, 64, 129]
    check v[63] and not v[62]
    v[63] = false
    v[64] = true      # already set: the count is unchanged
    check not v[63] and v.count == 3
    expect IndexDefect: discard v[130]

  test "logic and text forms":
    let a = parseBitVect("1100")
    let b = parseBitVect("1010")
    check $(a and b) == "1000"
    check $(a or b) == "1110"
    check (a and b).count == 1 and (a or b).count == 3
    check parseBitVect($a) == a
    expect ValueError: discard tanimoto(a, initBitVect(5))

  test "similarity":
    let a = parseBitVect("11110000")
    let b = parseBitVect("11001100")
    check tanimoto(a, b) == 2 / 6
    check dice(a, b) == 4 / 8
    check abs(cosine(a, b) - 0.5) < 1e-12
    check tversky(a, b, 1, 1) == tanimoto(a, b)
    check tversky(a, b, 0.5, 0.5) == dice(a, b)
    check tversky(parseBitVect("1100"), parseBitVect("1111"), 1, 0) == 1.0
    # empty fingerprints, as in RDKit
    check tanimoto(initBitVect(8), initBitVect(8)) == 0.0
    check dice(initBitVect(8), initBitVect(8)) == 0.0
    check bulkTanimoto(a, [a, b]) == @[1.0, 2 / 6]

suite "Morgan fingerprints":
  # bit IDs from RDKit (GetMorganGenerator(radius=2).GetSparseCountFingerprint)
  test "same IDs as RDKit":
    check morganCounts(parseSmiles("CCO")).ids == @[
      (864662311'u32, 1), (1535166686'u32, 1), (2245384272'u32, 1),
      (2246728737'u32, 1), (3542456614'u32, 1), (4018048386'u32, 1)]
    check morganCounts(parseSmiles("CCCO")).ids == @[
      (864662311'u32, 1), (1173125914'u32, 1), (1365892349'u32, 1),
      (1535166686'u32, 1), (2245384272'u32, 2), (2246728737'u32, 1),
      (3542456614'u32, 1), (4023654873'u32, 1)]
    check toSeq(morganFingerprint(parseSmiles("CCO"), nBits = 64).onBits) ==
          @[2, 16, 30, 33, 38, 39]

  test "radius 0 is the atom invariants":
    let m = parseSmiles("c1ccccc1O")
    var c = initCountTable[uint32]()
    for x in connectivityInvariants(m): c.inc x
    check morganCounts(m, radius = 0) == c

  test "similarity, compared with RDKit":
    let a = morganCounts(parseSmiles("CCO"))
    let b = morganCounts(parseSmiles("CCCO"))
    check tanimoto(a, b) == 0.5
    check abs(dice(a, b) - 2 / 3) < 1e-12
    let x = morganFingerprint(parseSmiles("CCO"), nBits = 64)
    let y = morganFingerprint(parseSmiles("CCCO"), nBits = 64)
    check abs(tversky(x, y, 0.3, 0.7) - 0.6756756756756758) < 1e-12

  test "feature invariants":
    # bits: donor 1, acceptor 2, aromatic 4, halogen 8, basic 16, acidic 32
    let m = parseSmiles("OC(=O)c1ccccc1")
    check featureInvariants(m) == @[1'u32, 32, 2, 4, 4, 4, 4, 4, 4]   # from RDKit
    check morganCounts(m, useFeatures = true) != morganCounts(m)

  test "independent of atom order and spelling":
    var r = initRand(3)
    for s in ["CC(=O)Oc1ccccc1C(=O)O", "C[C@H](N)C(=O)O", "c1ccc2[nH]ccc2c1", "[Na+].[Cl-]"]:
      let m = parseSmiles(s)
      for _ in 1..10:
        check morganCounts(m.shuffled(r), radius = 3) == morganCounts(m, radius = 3)
    check morganCounts(parseSmiles("C1=CC=CC=C1")) == morganCounts(parseSmiles("c1ccccc1"))
    check morganCounts(parseSmiles("[H]OC")) == morganCounts(parseSmiles("CO"))

  test "edge cases":
    check morganCounts(Mol()).len == 0
    check morganFingerprint(parseSmiles("[Na+].[Cl-]")).count == 2
    # isotopes count by whole-mass shift from the average mass, as in RDKit:
    # 13C (13.003 vs 12.011) doesn't change the invariant, 18O does
    check morganCounts(parseSmiles("[13CH4]")) == morganCounts(parseSmiles("C"))
    check morganCounts(parseSmiles("[18OH2]")).ids == @[(864666388'u32, 1)]
    expect ValueError: discard morganFingerprint(parseSmiles("C"), nBits = 0)

  test "isotopes and dummy atoms, as RDKit":
    # first-atom invariants from RDKit's GetConnectivityInvariants
    for (smiles, expected) in [("[1*]C(=O)O", 2342113506'u32), ("[16*]c1ccccc1", 2342113506'u32),
                               ("[51Cr]", 3829170480'u32), ("[177Lu]", 1931905113'u32),
                               ("[235U]", 529912953'u32), ("[99Tc]", 625568639'u32),
                               ("[32PH3]", 984478166'u32), ("[299C]", 2244242826'u32)]:
      checkpoint smiles
      check connectivityInvariants(parseSmiles(smiles))[0] == expected
    check abs(isotopeMass(24, 51) - 50.9447674) < 1e-9

  test "feature bits stop after 1,000 matching atoms, as in RDKit":
    let f = featureInvariants(parseSmiles(newSeqWith(200, "c1ccccc1").join(".")))
    check f.countIt((it and 4) != 0) == 1000
    check f[999] == 4 and f[1000] == 0

  test "usable from GC-safe code (threads)":
    proc work(s: string): int {.gcsafe.} =
      let m = parseSmiles(s)
      m.morganFingerprint.count + m.morganFingerprint(useFeatures = true).count +
        m.substructMatches(parseSmarts("[OH]")).len
    check work("c1ccccc1O") > 0
