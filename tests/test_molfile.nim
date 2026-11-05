import std/[unittest, strutils, sequtils, tables, os, random, math]
import chem
import ./helpers

const
  # RDKit's 2D and 3D molfiles for C[C@H](F)/C=C/Cl
  wedged2D =
    "\n" & """
     RDKit          2D

  6  5  0  0  0  0  0  0  0  0999 V2000
   -2.2500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -2.5981    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2500    1.2990    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  1
  2  3  1  0
  2  4  1  0
  4  5  2  0
  5  6  1  0
M  END
"""
  embedded3D =
    "\n" & """
     RDKit          3D

  6  5  0  0  0  0  0  0  0  0999 V2000
   -1.4131    0.4406   -0.1838 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4958   -0.6989    0.2418 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0261   -1.3844    1.3081 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.8518   -0.1387    0.5563 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6504    0.1830   -0.4417 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2318    0.8404   -0.0768 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0
  2  3  1  0
  2  4  1  0
  4  5  2  0
  5  6  1  0
M  END
"""

proc molText(
    atoms: openArray[string], bonds: openArray[string], props: openArray[string] = []
): string =
  ## A V2000 block with zero coordinates.
  result =
    "test\n  hand\n\n" & align($atoms.len, 3) & align($bonds.len, 3) &
    "  0  0  0  0  0  0  0  0999 V2000\n"
  for a in atoms:
    result &= "    0.0000    0.0000    0.0000 " & a & "\n"
  for b in bonds:
    result &= b & "\n"
  for p in props:
    result &= p & "\n"
  result &= "M  END\n"

proc errorLine(s: string): int =
  try:
    discard parseMolBlock(s)
    -1
  except MolfileError as e:
    e.line

suite "molfile reading":
  test "atoms, bonds and hydrogens":
    let m = parseMolBlock(
      molText(["C   0  0", "C   0  0", "O   0  0"], ["  1  2  1  0", "  2  3  1  0"])
    )
    check m.toSmiles == "CCO"
    check m.coords.len == 0 # all zero: no coordinates
    check m.molecularFormula == "C2H6O"

  test "charges, isotopes and atom maps":
    # atom-block charge codes (3 = +1, 5 = -1)
    check parseMolBlock(molText(["N   0  3", "O   0  5"], ["  1  2  1  0"])).toSmiles ==
      "[NH3+][O-]"
    # M  CHG replaces all atom-block charges
    check parseMolBlock(
      molText(["N   0  3", "O   0  5"], ["  1  2  1  0"], ["M  CHG  1   1   1"])
    ).toSmiles == "[NH3+]O"
    check parseMolBlock(
      molText(["C   0  0", "C   0  0"], ["  1  2  1  0"], ["M  ISO  1   2  13"])
    ).toSmiles == "C[13CH3]"
    check parseMolBlock(molText(["D   0  0", "C   0  0"], ["  1  2  1  0"])).toSmiles ==
      "[2H]C"
    let mapped = molText(["C   0  0  0  0  0  0  0  0  0  7  0  0"], [])
    check parseMolBlock(mapped).atoms[0].atomMap == 7

  test "explicit valence and radicals":
    # valence field 15 = no hydrogens
    check parseMolBlock(molText(["C   0  0  0  0  0 15"], [])).molecularFormula == "C"
    check parseMolBlock(molText(["C   0  0  0  0  0  2"], [])).molecularFormula == "CH2"
    check parseMolBlock(
      molText(["C   0  0", "C   0  0"], ["  1  2  1  0"], ["M  RAD  1   1   2"])
    ).molecularFormula == "C2H5"

  test "aromatic bonds":
    var bonds: seq[string]
    for i in 1 .. 6:
      bonds.add align($i, 3) & align($(i mod 6 + 1), 3) & "  4  0"
    check parseMolBlock(molText(newSeqWith(6, "C   0  0"), bonds)).toSmiles == "c1ccccc1"

  test "stereo from 2D wedges and 3D coordinates":
    for s in [wedged2D, embedded3D]:
      let m = parseMolBlock(s)
      check m.toSmiles == canonicalSmiles("C[C@H](F)/C=C/Cl")
      check m.coords.len == 6

  test "stereo from atom parities when there are no coordinates":
    let r = parseMolBlock(
      molText(
        ["C   0  0", "C   0  0  1", "F   0  0", "Cl  0  0"],
        ["  1  2  1  0", "  2  3  1  0", "  2  4  1  0"],
      )
    )
    let s = parseMolBlock(
      molText(
        ["C   0  0", "C   0  0  2", "F   0  0", "Cl  0  0"],
        ["  1  2  1  0", "  2  3  1  0", "  2  4  1  0"],
      )
    )
    check r.atoms[1].chiral != ctNone
    check r.toSmiles != s.toSmiles
    check r.toSmiles ==
      canonicalSmiles(
        s.toSmiles.replace("@@", "!").replace("@", "@@").replace("!", "@")
      )

  test "either and unmarked stereo":
    # bond stereo 3 on the double bond: no E/Z; no wedge: no chirality
    let s = wedged2D.replace("  4  5  2  0", "  4  5  2  3").replace(
        "  2  1  1  1", "  2  1  1  0"
      )
    check parseMolBlock(s).toSmiles == "CC(F)C=CCl"

  test "V3000":
    let s =
      "\n" & """
  test

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 C 0 0 0 0 -
M  V30 CHG=-1
M  V30 3 O 0 0 0 5 MASS=18
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(1 2)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
"""
    let m = parseMolBlock(s)
    check m.toSmiles == "C[CH-][18OH:5]"
    # VAL=0 is the default valence, VAL=-1 none
    check parseMolBlock(s.replace("CHG=-1", "VAL=0")).toSmiles == "CC[18OH:5]"
    check parseMolBlock(s.replace("CHG=-1", "VAL=-1")).toSmiles == "C[C][18OH:5]"

  test "property lines":
    # an atom value ("V  ") is one line: the CHG line after it still counts
    check parseMolBlock(
      molText(
        ["C   0  0", "O   0  0"],
        ["  1  2  1  0"],
        ["V    1 a value", "M  CHG  1   2  -1"],
      )
    ).toSmiles == "C[O-]"
    # an alias ("A  ") has its text on the next line
    check parseMolBlock(
      molText(
        ["C   0  0", "O   0  0"], ["  1  2  1  0"], ["A    1", "M  CHG  1   2  -1"]
      )
    ).toSmiles == "CO"
    # M  ISO replaces atom-block isotopes, but not D/T symbols
    check parseMolBlock(
      molText(["D   0  0", "C   0  0"], ["  1  2  1  0"], ["M  ISO  1   2  13"])
    ).toSmiles == "[2H][13CH3]"

  test "query atoms are rejected":
    check errorLine(molText(["A   0  0"], [])) == 5
    check errorLine(molText(["C   0  0"], [], ["M  ALS   1  2 F C   N   "])) == 6

  test "errors carry a line number":
    check errorLine(molText(["Xx  0  0"], [])) == 5
    check errorLine(molText(["C   0  0", "C   0  0"], ["  1  2  8  0"])) == 7
    check errorLine(molText(["C   0  0"], ["  1  2  1  0"])) == 6
    check errorLine(molText(["C   0  0"], []).replace("M  END\n", "")) > 0
    check errorLine("x\n") > 0

suite "molfile writing":
  test "round trips":
    for smi in [
      "CCO", "c1ccccc1O", "C[N+](C)(C)C.[Cl-]", "[13CH4]", "c1cc[nH]c1",
      "C[C@H](N)C(=O)O", "[CH3]", "[Na+].[O-]C(=O)C", "O=S1(=O)CCCC1",
      "C[S@](=O)c1ccccc1", "[2H]C([2H])([2H])O", "[CH2:3]=O",
    ]:
      let m = parseSmiles(smi)
      checkpoint smi
      check parseMolBlock(m.toMolBlock).toSmiles == m.toSmiles
      check parseMolBlock(m.toMolBlock(v3000 = true)).toSmiles == m.toSmiles
      check parseMolBlock(m.toMolBlock(kekule = false)).toSmiles == m.toSmiles

  test "without coordinates a layout keeps all stereo":
    let m = parseSmiles("C[C@H](F)/C=C/Cl")
    let back = parseMolBlock(m.toMolBlock)
    check back.toSmiles == m.toSmiles
    check back.coords.len == m.numAtoms and not back.is3D
    # without the layout: zero coordinates, E/Z lost, tetrahedral stereo
    # kept as parities
    check parseMolBlock(m.toMolBlock(layout = false)).toSmiles ==
      canonicalSmiles("C[C@H](F)C=CCl")

  test "with coordinates all stereo is kept":
    for s in [wedged2D, embedded3D]:
      let m = parseMolBlock(s)
      for v3 in [false, true]:
        let back = parseMolBlock(m.toMolBlock(v3000 = v3))
        check back.toSmiles == m.toSmiles
        check back.coords.len == 6
    # an unspecified double bond drawn in 2D is written as "either"
    let plain = parseMolBlock(
      wedged2D.replace("  2  1  1  1", "  2  1  1  0").replace(
        "  4  5  2  0", "  4  5  2  3"
      )
    )
    check "  4  5  2  3" in plain.toMolBlock
    check parseMolBlock(plain.toMolBlock).toSmiles == "CC(F)C=CCl"

  test "ring double bonds aren't marked as either":
    var benzene = parseSmiles("c1ccccc1")
    for k in 0 .. 5:
      let t = float(k) * 1.0471975512
      benzene.coords.add (cos(t), sin(t), 0.0)
    for line in benzene.toMolBlock.splitLines[10 .. 15]:
      check line.endsWith("  0")

  test "coordinates follow renumbering":
    var r = initRand(3)
    let m = parseMolBlock(wedged2D)
    let s = m.shuffled(r)
    check s.coords.len == 6
    check parseMolBlock(s.toMolBlock).toSmiles == m.toSmiles

  test "large molecules are written as V3000":
    let m = parseSmiles("C".repeat(1000))
    let s = m.toMolBlock
    check "V3000" in s
    check parseMolBlock(s).numAtoms == 1000

  test "header":
    let s = parseSmiles("C").toMolBlock("methane")
    check s.splitLines[0] == "methane"
    check "V2000" in s.splitLines[3]
    # dimension code in columns 21-22
    check parseMolBlock(wedged2D).toMolBlock.splitLines[1][20 .. 21] == "2D"
    check parseMolBlock(embedded3D).toMolBlock.splitLines[1][20 .. 21] == "3D"

  test "a chiral atom that can't be drawn raises":
    var m = parseSmiles("C[C@H](F)Cl")
    m.coords = newSeqWith(4, (1.0, 1.0, 0.0)) # all atoms on one point
    expect ChemError:
      discard m.toMolBlock
    m.coords = @[]
    check m.toMolBlock.len > 0 # without coordinates: parity

suite "SD files":
  const twoRecords =
    molText(["C   0  0"], []).replace("test", "methane") &
    "> <mw>\n16.04\n\n> <notes> (1)\nline one\nline two\n\n$$$$\n" &
    molText(["Xx  0  0"], []).replace("test", "bad") & "> <mw>\n0\n\n$$$$\n"

  test "records, names and data items":
    let recs = toSeq(sdfRecords(twoRecords))
    check recs.len == 2
    check recs[0].name == "methane"
    check recs[0].error == ""
    check recs[0].mol.toSmiles == "C"
    check toSeq(recs[0].props.keys) == @["mw", "notes"]
    check recs[0].props["notes"] == "line one\nline two"
    check recs[1].name == "bad"
    check "Xx" in recs[1].error
    check recs[1].props["mw"] == "0"

  test "writing and reading back":
    var text = ""
    for (smi, name) in [("CCO", "ethanol"), ("C[C@H](N)C(=O)O", "alanine")]:
      text &= parseSmiles(smi).toSdf(name, {"smiles": smi, "empty": ""})
    let recs = toSeq(sdfRecords(text))
    check recs.mapIt(it.name) == @["ethanol", "alanine"]
    check recs.mapIt(it.mol.toSmiles) == @["CCO", canonicalSmiles("C[C@H](N)C(=O)O")]
    check recs[1].props["smiles"] == "C[C@H](N)C(=O)O"
    check recs[1].props["empty"] == ""
    check toSeq(sdfRecords(recs[1].toSdf)).mapIt(it.props) == @[recs[1].props]

  test "data items the format can't hold":
    let m = parseSmiles("C")
    check toSeq(sdfRecords(m.toSdf(props = {"a": "x\ny\n"})))[0].props["a"] == "x\ny"
    for props in [@{"a": "x\n\ny"}, @{"a": "x\n$$$$"}, @{"a>b": "x"}, @{"a\nb": "x"}]:
      expect ValueError:
        discard m.toSdf(props = props)

  test "files, Windows line endings and a missing final $$$$":
    let path = getTempDir() / "chem_test.sdf"
    writeFile(path, twoRecords.replace("\n", "\r\n") & molText(["N   0  0"], []))
    defer:
      removeFile(path)
    let recs = toSeq(readSdf(path))
    check recs.len == 3
    check recs[0].props["notes"] == "line one\nline two"
    check recs[2].mol.toSmiles == "N"
