import std/[unittest, random, sequtils, strutils]
import chem
import ./helpers

proc hydrogenAtoms(m: Mol): int = m.atoms.countIt(it.element == 1)

suite "removing hydrogens":
  test "readers remove plain hydrogen atoms":
    check parseSmiles("[H]C([H])([H])[H]").toSmiles == "C"
    check parseSmiles("[H]O[H]").toSmiles == "O"
    check parseSmiles("[H]C(=O)O").toSmiles == "O=CO"
    check parseSmiles("[H][N+](C)(C)C").toSmiles == "C[NH+](C)C"
    check parseSmiles("[H][Fe]").toSmiles == "[FeH]"
    check parseSmiles("[H]C", removeHs = false).hydrogenAtoms == 1
    check parseSmiles("[H]C", sanitize = false).hydrogenAtoms == 1

  test "H counts that the valence model wouldn't give become explicit":
    let m = parseSmiles("[H]n1cccc1")
    check m.toSmiles == "c1cc[nH]c1"
    check m.hydrogenAtoms == 0

  test "hydrogens that carry information are kept":
    for s in ["[2H]C", "[H][H]", "[H]", "[H+]", "[H-]", "*[H]", "[H:2][CH3:1]",
              "[H]/N=C/C", "[H]/C=C/F"]:
      checkpoint s
      check parseSmiles(s).hydrogenAtoms > 0
    # a bridging H (fails the valence check, so read unsanitized)
    check parseSmiles("B1[H]B[H]1", sanitize = false).removeHs.hydrogenAtoms == 2

  test "chirality is re-expressed":
    check parseSmiles("[H][C@](F)(Cl)Br").toSmiles == "F[C@@H](Cl)Br"
    check parseSmiles("F[C@]([H])(Cl)Br").toSmiles == "F[C@H](Cl)Br"
    check parseSmiles("F[C@@](Cl)(Br)[H]").toSmiles ==
          parseSmiles("F[C@@H](Cl)Br").toSmiles
    check parseSmiles("[H][C@@]1(C)CC[C@@H](C)CC1").toSmiles == "C[C@H]1CC[C@H](C)CC1"
    # a three-coordinate centre would lose a position: the H stays
    check parseSmiles("C[S@]([H])=O").toSmiles == "[H][S@@](C)=O"

  test "double-bond stereo moves to the other substituent":
    check parseSmiles("[H]/C(Cl)=C/F").toSmiles == parseSmiles("Cl/C=C\\F").toSmiles
    check parseSmiles("F/C([H])=C(/[H])Cl").toSmiles == parseSmiles("F/C=C\\Cl").toSmiles

  test "unsanitized molecules":
    var m = parseSmiles("[H]C([H])=O", sanitize = false)
    check m.removeHs.hydrogenAtoms == 0
    check m.removeHs.numAtoms == 2

suite "adding hydrogens":
  test "every H becomes an atom":
    let m = parseSmiles("CC(=O)Nc1cc[nH]c1").addHs
    check m.hydrogenAtoms == 8
    check m.atoms.allIt(it.explicitH == 0)
    check m.molecularFormula == parseSmiles("CC(=O)Nc1cc[nH]c1").molecularFormula
    check (0..<m.numAtoms).toSeq.allIt(m.atoms[it].element == 1 or totalHCount(m, it) == 0)

  test "round trips keep the molecule, stereo included":
    var r = initRand(7)
    for s in ["C[C@H](F)Cl", "F/C=C/F", "C[C@H]1CC[C@@H](C)CC1", "c1cc[nH]c1",
              "C[S@](=O)c1ccccc1", "[NH4+]", "[CH2]C", "O[C@H]1[C@H](O)[C@@H](O)CC1",
              "N[C@@H](Cc1ccccc1)C(=O)O", "C/C=C/C=C\\C", "[13CH3][C@H](O)F"]:
      checkpoint s
      let m = parseSmiles(s)
      let h = m.addHs
      check h.removeHs.toSmiles == m.toSmiles
      check parseSmiles(h.toSmiles).toSmiles == m.toSmiles
      check parseSmiles(h.shuffled(r).toSmiles(canonical = false)).toSmiles == m.toSmiles

  test "coordinates for new atoms, 2D":
    for text in ["C[C@H](F)CC", "C[C@@H](O)C(C)C", "N[C@@](C)(F)CO"]:
      checkpoint text
      var m = parseSmiles(text)
      # a zigzag chain, with atom 2 (a branch) below atom 1
      m.coords = (0..<m.numAtoms).toSeq.mapIt((it.float * 1.3, (it mod 2).float * 0.75, 0.0))
      m.coords[2] = (1.3, -1.5, 0.0)
      let h = m.addHs
      check h.coords.len == h.numAtoms
      check (m.numAtoms..<h.numAtoms).toSeq.allIt(h.coords[it].z == 0.0)
      # wedges drawn to explicit H atoms read back as the same stereo
      check parseMolBlock(h.toMolBlock).toSmiles == m.toSmiles

  test "coordinates for new atoms, 3D":
    var m = parseSmiles("FC(Cl)Br")
    m.coords = @[(1.0, 1.0, 1.0), (0.0, 0.0, 0.0), (1.0, -1.0, -1.0), (-1.0, 1.0, -1.0)]
    m.stereoFromCoords([bwNone, bwNone, bwNone], [false, false, false])
    check m.atoms[1].chiral != ctNone
    let h = m.addHs
    let p = h.coords[4]
    check p.x < 0 and p.y < 0 and p.z > 0   # the free corner
    check parseMolBlock(h.toMolBlock).toSmiles == m.toSmiles

  test "molfiles with hydrogen atoms":
    let m = parseSmiles("N[C@@H](C)C(=O)O")
    let text = m.addHs.toMolBlock
    check text.splitLines[3].strip.startsWith("13")   # 6 heavy atoms, 7 H
    check parseMolBlock(text).toSmiles == m.toSmiles
    check parseMolBlock(text, removeHs = false).hydrogenAtoms == 7
