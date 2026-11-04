import std/[unittest, random, algorithm, sequtils, strutils]
import chem
import ./helpers

proc count(smiles, smarts: string; uniquify = true): int =
  parseSmiles(smiles).substructMatches(parseSmarts(smarts), uniquify = uniquify).len

proc hits(smiles, smarts: string): bool =
  parseSmiles(smiles).hasSubstructMatch(parseSmarts(smarts))

proc chiralHit(smiles, query: string): bool =
  parseSmiles(smiles).hasSubstructMatch(parseSmarts(query), useChirality = true)

suite "SMARTS parsing":
  test "errors carry a position":
    for (s, pos) in [("C(C", 3), ("[C", 0), ("C1CC", 4), ("CC=", 2), ("[Zz]", 1),
                     ("C)", 1), ("[$(C]", 1), ("", 0), ("[C;$(C=[O;Q])]", 10),
                     ("C/C=C(/F)/Cl", 9), ("C-1CC=1", 6), ("[C@SP1](F)(Cl)Br", 2)]:
      try:
        discard parseSmarts(s)
        checkpoint s
        check false
      except SmartsError as e:
        checkpoint s & ": " & e.msg
        check e.pos == pos

  test "atom maps are kept":
    let q = parseSmarts("[C:1](=O)[O;H1:2]")
    check q.atoms.mapIt(it.atomMap) == @[1, 0, 2]

  test "patterns from the RDKit data files parse":
    for s in ["[CH1X4]([N,O,P,S,F,Cl,Br,I])([A;!#1])[A;!#1]",
              "[NH0;+,+2,+3](=[A;!#1])([A;!#1])[!#1;A,a]",
              "[$(a-[CH3]),$(a-[CH2]-[CH3]),$(a-[CH2]-[CH2]~[!N;!O]);!$(a(:a!:*):a!:*)]",
              "c:1:c:c(:c:c:c:1-[#6;X4]-c:2:c:c:c(:c:c:2)-[#7&H2,$([#7;!H0]-[#6;X4])])",
              "[#6]-[#16;X2]-[#6](=[#16])", "[OX1;-,-2,-3][!#1;!N;!S]", "C-,:C", "*!@*"]:
      checkpoint s
      check parseSmarts(s).numAtoms > 0

  test "chirality followed by an uppercase primitive":
    check chiralHit("C[C@H](F)Cl", "[C@D3](C)(F)Cl")
    check not chiralHit("C[C@H](F)Cl", "[C@@D3](C)(F)Cl")
    check chiralHit("C[C@H](F)Cl", "[C@X4](C)(F)Cl")
    check chiralHit("C[C@H](F)Cl", "[C@TH1H](C)(F)Cl") == chiralHit("C[C@H](F)Cl", "[C@H](C)(F)Cl")

  test "ring closure bonds":
    check count("C1CC1", "C-1CC-1") == 1
    check count("C1CC1", "C-1CC1") == 1
    check count("C/C1=C/CC1", "C1CC/1") == count("C/C1=C/CC1", "C1CC1")

suite "atom primitives":
  test "elements and aromaticity":
    check count("c1ccccc1C", "c") == 6
    check count("c1ccccc1C", "C") == 1
    check count("c1ccccc1C", "[#6]") == 7
    check count("c1ccccc1C", "a") == 6
    check count("c1ccccc1C", "A") == 1
    check count("c1ccncc1", "[n]") == 1
    check count("C[Se]C.c1cc[se]c1", "[se]") == 1
    check count("C[Se]C.c1cc[se]c1", "[Se]") == 1
    check count("ClCCBr", "[Cl,Br]") == 2
    check count("CC", "*") == 2

  test "hydrogens and connectivity":
    let s = "CC(C)(C)CO"                 # neopentyl alcohol
    check count(s, "[CH3]") == 3
    check count(s, "[CH2]") == 1
    check count(s, "[CH0]") == 1
    check count(s, "[D4]") == 1
    check count(s, "[D1]") == 4
    check count(s, "[X4]") == 5
    check count(s, "[OX2H]") == 1
    check count("C[NH3+]", "[NX4+]") == 1
    check count("C[nH]1cccc1", "[nH]") == 0   # not a valid pyrrole spelling
    check count("c1cc[nH]c1", "[nH]") == 1
    check count("c1cc[nH]c1", "[h]") == 5      # any H not stored as an atom
    # hydrogen atoms are removed when reading, unless asked not to
    check count("[H]C([H])C", "[#1]") == 0
    check count("[H]C([H])C", "[H3]") == 2
    let withHs = parseSmiles("[H]C([H])C", removeHs = false)
    for (smarts, n) in [("[h1]", 1), ("[H3]", 2), ("[#1]", 2), ("[H]", 2)]:
      check withHs.substructMatches(parseSmarts(smarts)).len == n
    check count("[2H]C", "[2H]") == 1
    check count("[2H]C", "[2#1]") == 1
    check count("[H+].[Cl-]", "[H+]") == 1

  test "valence, charge, isotope":
    check count("c1ccccc1O", "[v4]") == 6
    check count("CS(=O)(=O)C", "[Sv6]") == 1
    check count("C[N+](C)(C)C.[O-]C", "[+]") == 1
    check count("C[N+](C)(C)C.[O-]C", "[-1]") == 1
    check count("C[N+](C)(C)C.[O-]C", "[+0]") == 5
    check count("[13CH4].C", "[13C]") == 1
    check count("[Fe+3]", "[Fe+++]") == 1

  test "rings":
    let s = "c1ccc2ccccc2c1CC1CC1"
    check count(s, "[R]") == 13
    check count(s, "[R0]") == 1
    check count(s, "[R2]") == 2
    check count(s, "[r3]") == 3
    check count(s, "[r6]") == 10
    check count(s, "[x3]") == 2
    check count(s, "[x]") == 13
    check count(s, "[!R]") == 1
    check count("C1CC2CCC1CC2", "[R3]") == 2   # bicyclo[2.2.2]octane bridgeheads

  test "logical operators and precedence":
    let s = "CNC(=O)O"
    check count(s, "[C,N;H1]") == 1       # (C or N) and H1
    check count(s, "[C,N&H1]") == 3       # C or (N and H1)
    check count(s, "[!C]") == 3
    check count(s, "[!C;!N]") == 2
    check count(s, "[C&!H3]") == 1
    check count(s, "[!!C]") == 2

  test "recursive SMARTS":
    let s = "CC(=O)OCC(=O)O"
    check count(s, "[$(C=O)]") == 2
    check count(s, "[$(C(=O)[OH])]") == 1
    check count(s, "[O;$(OC=O);!$(O[H]);!$([OH])]") == 1
    check count(s, "[C;!$(C=O)]") == 2
    check count(s, "[$([OX2]C(=O)[$([CH3])])]") == 1      # nested
    check count(s, "[$(*=O)]~[O;!H1]") == 3

suite "bonds and topology":
  test "bond primitives":
    check count("c1ccccc1C=CC#N", "c:c") == 6
    check count("c1ccccc1C=CC#N", "c-c") == 0
    check count("c1ccccc1C=CC#N", "cC") == 1      # unwritten: single or aromatic
    check count("c1ccccc1C=CC#N", "cc") == 6
    check count("c1ccccc1C=CC#N", "C=C") == 1
    check count("c1ccccc1C=CC#N", "C#N") == 1
    check count("c1ccccc1C=CC#N", "*~*") == 10
    check count("c1ccccc1C=CC#N", "*@*") == 6
    check count("c1ccccc1C=CC#N", "*!@*") == 4
    check count("c1ccccc1C=CC#N", "[#6]=,#[#6,#7]") == 2
    check count("C1CCCCC1C", "C@C") == 6
    check count("C1CCCCC1C", "C-;!@C") == 1

  test "ring closures, branches and disconnected queries":
    check count("C1CCCCC1", "C1CCCCC1") == 1
    check count("C1CCCC1", "C1CCCCC1") == 0
    check count("OCC(=O)O", "C(=O)O") == 1
    check count("CCO.CCN", "O.N") == 1
    check count("CCO", "O.N") == 0
    check count("CC", "C.C") == 1

  test "uniquify and maxMatches":
    check count("c1ccccc1", "c1ccccc1") == 1
    check count("c1ccccc1", "c1ccccc1", uniquify = false) == 12
    check count("CC(C)(C)C", "C(C)(C)C", uniquify = false) == 24   # 4!
    check count("CC(C)(C)C", "C(C)(C)C") == 4     # distinct atom sets
    let m = parseSmiles("CCCCCCCC")
    check m.substructMatches(parseSmarts("C"), maxMatches = 3).len == 3

  test "matches map query atoms to molecule atoms":
    let m = parseSmiles("OC(=O)c1ccccc1")
    let match = m.substructMatch(parseSmarts("c-C(=O)[OH]"))
    check match.len == 4
    check m.atoms[match[0]].aromatic
    check m.atoms[match[3]].element == 8
    check m.bondBetween(match[0], match[1]) >= 0
    check m.substructMatch(parseSmarts("N")).len == 0

  test "results don't depend on atom order":
    var r = initRand(5)
    let m = parseSmiles("CC(=O)Nc1ccc(O)cc1C1CC1")
    for smarts in ["[$(C=O)]N", "c[OH]", "[R]", "[r3]", "*", "C(=O)N", "a:a"]:
      let q = parseSmarts(smarts)
      let n = m.substructMatches(q).len
      for _ in 0..<10:
        check m.shuffled(r).substructMatches(q).len == n

suite "molecule queries":
  test "SMILES as a query":
    let aspirin = parseSmiles("CC(=O)Oc1ccccc1C(=O)O")
    check aspirin.hasSubstructMatch(parseSmiles("c1ccccc1"))
    check aspirin.hasSubstructMatch(parseSmiles("CC=O"))      # H counts not compared
    check aspirin.hasSubstructMatch(parseSmiles("C1=CC=CC=C1C"))  # Kekulé input is aromatic again
    check not aspirin.hasSubstructMatch(parseSmiles("C1=CCCC=C1"))
    check aspirin.substructMatches(parseSmiles("C(=O)O")).len == 2
    check not parseSmiles("CC(=O)[O-]").hasSubstructMatch(parseSmiles("C[O+]"))
    check parseSmiles("CC(=O)[O-]").hasSubstructMatch(parseSmiles("[O-]"))
    check parseSmiles("CC(=O)[O-]").hasSubstructMatch(parseSmiles("O"))   # charge 0: any

  test "string convenience overloads":
    let m = parseSmiles("CCO")
    check m.hasSubstructMatch("[OH]")
    check m.substructMatches("C").len == 2

suite "stereo-aware matching":
  test "tetrahedral centres":
    let r = "C[C@@H](N)C(=O)O"     # L-alanine
    let s = "C[C@H](N)C(=O)O"
    check chiralHit(r, "C[C@@H](N)C(=O)O")
    check not chiralHit(r, "C[C@H](N)C(=O)O")
    check chiralHit(r, "N[C@H](C)C(=O)O")     # same centre, other order
    check chiralHit(r, "C[C@@H](N)C")         # partial
    check chiralHit(r, "C[C@@](N)C=O")        # H not written: it goes last
    check not chiralHit(r, "C[C@](N)C=O")
    check chiralHit(r, "C[C@](N)C")           # either C can be either C
    check chiralHit(s, "C[C@H](N)C")
    check hits(r, "C[C@H](N)C(=O)O")          # ignored without useChirality
    check not chiralHit("CC(N)C(=O)O", "C[C@@H](N)C(=O)O")   # unspecified target
    check chiralHit("CC(N)C(=O)O", "CC(N)C(=O)O")
    check chiralHit("[H][C@](C)(N)C(=O)O", "[C@@H](C)(N)C(=O)O") ==
          chiralHit(canonicalSmiles("[H][C@](C)(N)C(=O)O"), "[C@@H](C)(N)C(=O)O")

  test "molecule queries with stereo":
    let r = parseSmiles("C[C@@H](N)C(=O)O")
    check r.hasSubstructMatch(parseSmiles("N[C@H](C)C(=O)O"), useChirality = true)
    check not r.hasSubstructMatch(parseSmiles("N[C@@H](C)C(=O)O"), useChirality = true)

  test "double bonds":
    let e = "C/C=C/C(=O)O"      # (E)-crotonic acid
    check chiralHit(e, "C/C=C/C")
    check chiralHit(e, "C\\C=C\\C")
    check chiralHit(e, "C(\\C)=C/C")
    check not chiralHit(e, "C(/C)=C/C")
    check not chiralHit(e, "C/C=C\\C")
    check chiralHit("C/C=C\\C(=O)O", "C/C=C\\C")
    check hits(e, "C/C=C\\C")              # ignored without useChirality
    check not chiralHit("CC=CC(=O)O", "C/C=C/C")
    check chiralHit("F/C(Cl)=C/F", "Cl/C=C\\F")   # via the other substituent

  test "stereo inside recursive SMARTS":
    check chiralHit("C[C@H](F)Cl", "[$(C[C@H](F)Cl)]")
    check not chiralHit("C[C@H](F)Cl", "[$(C[C@@H](F)Cl)]")
    check hits("C[C@H](F)Cl", "[$(C[C@@H](F)Cl)]")        # ignored without useChirality
    # read as outside $(): the H of a first atom comes first (RDKit puts it
    # second inside recursive SMARTS)
    check chiralHit("C[C@H](F)Cl", "[$([C@@H](C)(F)Cl)]") ==
          chiralHit("C[C@H](F)Cl", "[C@@H](C)(F)Cl")
    check chiralHit("C/C=C/C", "[$(C/C=C/C)]")
    check not chiralHit("C/C=C/C", "[$(C/C=C\\C)]")
