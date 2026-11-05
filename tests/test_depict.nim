import std/[unittest, math, random, sequtils]
import chem
import ./helpers

func dist(m: Mol, i, j: int): float =
  hypot(m.coords[i].x - m.coords[j].x, m.coords[i].y - m.coords[j].y)

proc bondRange(m: Mol): (float, float) =
  result = (Inf, 0.0)
  for bd in m.bonds:
    let d = dist(m, bd.a, bd.b)
    result = (min(result[0], d), max(result[1], d))

proc closestNonBonded(m: Mol): float =
  result = Inf
  for i in 0 ..< m.numAtoms:
    for j in i + 1 ..< m.numAtoms:
      if m.bondBetween(i, j) < 0:
        result = min(result, dist(m, i, j))

proc sameSide(m: Mol, bi: int): bool =
  ## whether a stereo double bond's reference atoms are drawn on one side
  let bd = m.bonds[bi]
  let (a, b) = (m.coords[bd.a], m.coords[bd.b])
  proc side(p: Point3): float =
    (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x)

  side(m.coords[bd.stereoAtoms[0]]) * side(m.coords[bd.stereoAtoms[1]]) > 0

suite "2D layout":
  test "regular rings and chains":
    for s in [
      "c1ccccc1", "c1ccc2ccccc2c1", "c1ccc2cccc2cc1", "c1ccccc1-c1ccccc1", "CCCCCCCCCC",
      "CC(C)(C)C(=O)NC", "C[C@H](N)C(=O)O", "C1CCC2(CC1)CCCC2",
      "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@@]34C)[C@@H]1CC[C@@H]2O",
      "C1COCCOCCOCCOCCOCCO1", "C#CC#CC", "C=C=C", "O=C(O)c1ccccc1OC(C)=O",
      "c1ccc2c(c1)c1ccccc1c1ccccc21", "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    ]:
      checkpoint s
      let m = parseSmiles(s).with2DCoords
      check m.coords.len == m.numAtoms and m.coords.allIt(it.z == 0)
      let (lo, hi) = bondRange(m)
      check abs(lo - bondLength) < 1e-6 and abs(hi - bondLength) < 1e-6
      check closestNonBonded(m) > 0.4 * bondLength

  test "cages and bridged rings keep atoms apart":
    for s in ["C1CC2CCC1C2", "C1C2CC3CC1CC(C2)C3", "C12C3C4C1C5C2C3C45"]:
      checkpoint s
      let m = parseSmiles(s).with2DCoords
      let (lo, hi) = bondRange(m)
      check lo > 0.6 * bondLength and hi < 1.6 * bondLength
      check closestNonBonded(m) > 0.3 * bondLength

  test "triple bonds and allenes are straight":
    let m = parseSmiles("CC#CC").with2DCoords
    check abs(dist(m, 0, 3) - 3 * bondLength) < 1e-6

  test "the same drawing for any atom order":
    var r = initRand(7)
    for s in [
      "CC(=O)Oc1ccccc1C(=O)O", "C/C=C/C(Cl)=C\\Br", "C[C@H]1CC[C@@H](C)CC1",
      "C1CC2CCC1C2.[Na+]", "c1ccc2[nH]ccc2c1",
    ]:
      checkpoint s
      let m = parseSmiles(s)
      let a = m.with2DCoords
      let ranks = canonicalRanks(m)
      for _ in 1 .. 5:
        let shuffledMol = m.shuffled(r)
        let b = shuffledMol.with2DCoords
        let rb = canonicalRanks(shuffledMol)
        for i in 0 ..< m.numAtoms:
          let j = rb.find(ranks[i])
          check abs(a.coords[i].x - b.coords[j].x) < 1e-9
          check abs(a.coords[i].y - b.coords[j].y) < 1e-9

  test "cis/trans as stored, also in large rings":
    for s in [
      "F/C=C/F", "F/C=C\\F", "C/C=C/C=C\\C", "CC/C(C)=C(/CC)C(C)C", "C1/C=C\\CCCCC1",
      "C1/C=C/CCCCC1", "C1CCCCC/C=C/C=C/CCCCC1", "C1CCCC/C=C\\C=C/C=C\\CCCCC1",
      "OC1=CC/C(=C(/C)C2CCCC2)C1",
    ]:
      checkpoint s
      let m = parseSmiles(s).with2DCoords
      for bi, bd in m.bonds:
        if bd.stereo != bsNone:
          check sameSide(m, bi) == (bd.stereo == bsCis)

  test "stereo survives drawing and reading back":
    for s in [
      "C[C@H](N)C(=O)O", "C[C@H](F)/C=C/Cl", "C[C@H]1CC[C@@H](C)CC1",
      "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "C/C=C\\[C@@H](C)CC",
      "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@@]34C)[C@@H]1CC[C@@H]2O",
    ]:
      checkpoint s
      let m = parseSmiles(s)
      check parseMolBlock(m.toMolBlock).toSmiles == m.toSmiles
      check parseMolBlock(m.with2DCoords.toMolBlock(v3000 = true)).toSmiles == m.toSmiles

  test "fragments side by side, edge cases":
    let m = parseSmiles("[Na+].[Cl-].CCO").with2DCoords
    check m.coords[0].x < m.coords[1].x and m.coords[1].x < m.coords[2].x
    check Mol().with2DCoords.coords.len == 0
    check parseSmiles("C").with2DCoords.coords == @[(0.0, 0.0, 0.0)]
    # replaces 3D coordinates
    var t = parseSmiles("CO")
    t.coords = @[(0.0, 0.0, 1.0), (1.0, 1.0, 1.0)]
    t.compute2DCoords
    check t.coords.allIt(it.z == 0) and abs(dist(t, 0, 1) - bondLength) < 1e-9
