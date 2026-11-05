import std/[unittest, random, sequtils, strutils]
import chem
import ./helpers

# expected strings from the reference `selfies` package (2.x)

suite "SELFIES encoding":
  test "same as the reference encoder for Kekulé SMILES":
    for (smiles, expected) in [
        ("C=CF", "[C][=C][F]"),
        ("C1=CC=CC=C1", "[C][=C][C][=C][C][=C][Ring1][=Branch1]"),
        ("CC(=O)OC1=CC=CC=C1C(=O)O",
         "[C][C][=Branch1][C][=O][O][C][=C][C][=C][C][=C][Ring1][=Branch1][C][=Branch1][C][=O][O]"),
        ("C[C@@H](N)C(=O)O", "[C][C@@H1][Branch1][C][N][C][=Branch1][C][=O][O]"),
        ("F/C=C/F", "[F][/C][=C][/F]"),
        ("C1CC2CC1C2", "[C][C][C][C][C][Ring1][Branch1][C][Ring1][Ring2]"),
        ("[NH4+].[Cl-]", "[NH4+1].[Cl-1]"),
        ("[13CH3]O", "[13CH3][O]"),
        ("C1CC[C@H]2CCCC[C@@H]2C1",
         "[C][C][C][C@H1][C][C][C][C][C@@H1][Ring1][=Branch1][C][Ring1][#Branch2]"),
        ("N#CC=C", "[N][#C][C][=C]"),
        ("O=S(=O)(O)O", "[O][=S][=Branch1][C][=O][Branch1][C][O][O]")]:
      checkpoint smiles
      check smilesToSelfies(smiles) == expected

  test "aromatic SMILES is kekulized first":
    check smilesToSelfies("c1ccccc1") == "[C][=C][C][=C][C][=C][Ring1][=Branch1]"
    check smilesToSelfies("CC[Co]") == "[C][C][Co]"   # not aromatic

  test "canonical: the same for any atom order":
    var r = initRand(11)
    for s in ["CC(=O)Oc1ccccc1C(=O)O", "C[C@H](N)C(=O)O", "F/C=C/C(Cl)=C\\Br",
              "C[C@H]1CC[C@@H](C)CC1", "c1ccc2[nH]ccc2c1", "[Na+].[O-]C(=O)C"]:
      checkpoint s
      let m = parseSmiles(s)
      for _ in 1..10:
        check m.shuffled(r).toSelfies == m.toSelfies

  test "molecules SELFIES can't hold":
    expect SelfiesError: discard parseSmiles("ClI(Cl)Cl").toSelfies   # I has 3 bonds
    check parseSmiles("ClI(Cl)Cl").toSelfies(constraints = scHypervalent).len > 0
    check parseSmiles("ClI(Cl)Cl").toSelfies(strict = false).len > 0
    expect SelfiesError: discard parseSmiles("C*").toSelfies

  test "invalid SMILES":
    for smiles in ["", "C11", "[C:]O", "[99999999999999999999CH4]", "C(", "C1CC"]:
      checkpoint smiles
      expect SelfiesError: discard smilesToSelfies(smiles)

  test "aromatic input keeps its atoms in order":
    # from the reference encoder
    check smilesToSelfies("[H]c1ccccc1") == "[H][C][=C][C][=C][C][=C][Ring1][=Branch1]"

  test "deeply nested branches":
    proc nested(n: int): string = "C(".repeat(n) & "C" & ")C".repeat(n)
    let sf = smilesToSelfies(nested(500))
    check parseSelfies(sf).numAtoms == 1001
    # a branch holds at most 4,095 symbols, so nesting is bounded
    expect SelfiesError: discard smilesToSelfies(nested(5000))

suite "SELFIES decoding":
  test "same SMILES as the reference decoder":
    for (sf, expected) in [
        ("[C][=C][F]", "C=CF"),
        ("[C][Branch1][C][O][O]", "C(O)O"),
        ("[C][O][=O]", "COO"),                # bond orders lowered to fit
        ("[F][=C][=C][#N]", "FC=C=N"),
        ("[C][C][C][Ring1][Ring1]", "C1CC1"),
        ("[C][C][C][=Ring1][Ring1]", "C=1CC=1"),
        ("[Branch1][C][C]", "CC"),            # nothing to branch from
        ("[C][epsilon][O]", "C"),
        ("[C][nop][O]", "CO"),
        ("[O][=O][C]", "O=O"),                # O is full
        ("[N][Branch1][C][C][Branch1][C][C][Branch1][C][C][C]", "N(C)(C)CCC"),
        ("[CH4][C]", "[CH4]"),
        ("[C].[O]", "C.O"),
        ("[S][=C][=Ring1][Ring1]", "S#C"),    # a ring bond onto a bond
        ("[C][/C][=C][/C]", "C/C=C/C"),
        ("[C][C@@H1][Branch1][C][O][N]", "C[C@@H1](O)N")]:
      checkpoint sf
      check selfiesToSmiles(sf) == expected
    check selfiesToSmiles("[S][=O][=O][=O]", scOctetRule) == "S=O"
    check selfiesToSmiles("[I][Branch1][C][F][Branch1][C][F][F]", scHypervalent) == "I(F)(F)F"

  test "invalid symbols":
    for sf in ["[Xx]", "[C][Branch4]", "[CH]", "[C+0]", "[CH9]", "[C",
               "[99999999999999999999C]", "[C][C+99999999999999999999]"]:
      checkpoint sf
      expect SelfiesError: discard selfiesToSmiles(sf)

  test "round trips":
    for s in ["CC(=O)Oc1ccccc1C(=O)O", "C[C@H](N)C(=O)O", "F/C=C/C(Cl)=C\\Br",
              "C[C@H]1CC[C@@H](C)CC1", "c1ccc2[nH]ccc2c1", "[Na+].[O-]C(=O)C",
              "C[S@](=O)c1ccccc1", "[2H]C([2H])([2H])O", "C1CC2(CC1)CCCC2",
              "C12C3C4C1C5C2C3C45"]:
      checkpoint s
      let m = parseSmiles(s)
      check parseSelfies(m.toSelfies).toSmiles == m.toSmiles
      check parseSelfies(m.toSelfies(canonical = false)).toSmiles == m.toSmiles

  test "every string of robust symbols decodes to a molecule":
    let alphabet = robustAlphabet()
    check "[=N+1]" in alphabet and "[Branch3]" in alphabet and "[#O]" notin alphabet
    var r = initRand(5)
    for _ in 1..500:
      let s = newSeqWith(r.rand(1..40), r.sample(alphabet)).join()
      checkpoint s
      discard parseSelfies(s)

  test "deeply nested branches":
    # every branch gets its own budget, so nesting is unbounded
    let s = "[C][Branch3][P][P][P]".repeat(200000)
    check selfiesToSmiles(s).len > 0

  test "symbols":
    check toSeq(selfiesSymbols("[C][=C][F].[C]")) == @["[C]", "[=C]", "[F]", ".", "[C]"]
