import std/[unittest, algorithm]
import chem

proc ringSizes(s: string): seq[int] =
  for r in parseSmiles(s).sssr: result.add r.len
  result.sort()

proc aromaticAtoms(s: string): int =
  for a in parseSmiles(s).atoms:
    if a.aromatic: inc result

suite "SSSR":
  test "ring sizes":
    check ringSizes("CCO").len == 0
    check ringSizes("c1ccccc1") == @[6]
    check ringSizes("c1ccc2ccccc2c1") == @[6, 6]
    check ringSizes("C12C3C4C1C5C2C3C45") == @[4, 4, 4, 4, 4]   # cubane
    check ringSizes("C1CC2CCC1CC2") == @[6, 6]                  # bicyclo[2.2.2]octane
    check ringSizes("C1CCC2(CC1)CCCC2") == @[5, 6]              # spiro
    check ringSizes("C1CC1.C1CCC1") == @[3, 4]

  test "symmetrized SSSR":
    proc symmSizes(s: string): seq[int] =
      for r in parseSmiles(s).symmetrizedSssr: result.add r.len
      result.sort()
    check symmSizes("CCO").len == 0
    check symmSizes("c1ccc2ccccc2c1") == @[6, 6]
    check symmSizes("C12C3C4C1C5C2C3C45") == @[4, 4, 4, 4, 4, 4]   # cubane
    check symmSizes("C1CC2CCC1CC2") == @[6, 6, 6]
    check symmSizes("C1C2CC3CC1CC(C2)C3") == @[6, 6, 6, 6]         # adamantane
    check symmSizes("C1CC2CCC1C2") == @[5, 5]                      # norbornane
    check symmSizes("C1CC1.C1CCC1") == @[3, 4]

  test "rings are in walking order":
    let m = parseSmiles("c1ccc2ccccc2c1")
    for ring in m.sssr:
      for k in 0..<ring.len:
        check m.bondBetween(ring[k], ring[(k + 1) mod ring.len]) >= 0

suite "kekulization":
  test "kekulize benzene":
    var m = parseSmiles("c1ccccc1")
    kekulize(m)
    var doubles = 0
    for bd in m.bonds:
      check bd.order in {boSingle, boDouble}
      if bd.order == boDouble: inc doubles
    check doubles == 3
    for i in 0..<m.numAtoms:
      check not m.atoms[i].aromatic
      check totalHCount(m, i) == 1

  test "kekule SMILES output":
    check parseSmiles("c1ccccc1").toSmiles(kekule = true) == "C1=CC=CC=C1"
    check parseSmiles("c1ccncc1").toSmiles(kekule = true) == "C1=CC=NC=C1"
    # odd-ring fused system needs a general (blossom) matching
    check canonicalSmiles(parseSmiles("c1ccc2cccc2cc1").toSmiles(kekule = true)) ==
          canonicalSmiles("c1ccc2cccc2cc1")

  test "kekule output is canonical":
    let a = parseSmiles("c1ccc2ccccc2c1").toSmiles(kekule = true)
    let b = parseSmiles("C1=CC2=CC=CC=C2C=C1").toSmiles(kekule = true)
    check a == b

  test "impossible structures raise KekulizeError":
    for bad in ["c1cccc1", "c1ccnc1", "c1ccccc1(C)C", "c1cc(=O)ccc1"]:
      expect KekulizeError:
        discard parseSmiles(bad)

suite "aromaticity perception":
  test "Kekule and aromatic input agree":
    for (kek, arom) in [("C1=CC=CC=C1", "c1ccccc1"),
                        ("C1=CC=NC=C1", "c1ccncc1"),
                        ("C1=CC=CN1", "c1cc[nH]c1"),
                        ("C1=COC=C1", "c1ccoc1"),
                        ("C1=CSC=C1", "c1ccsc1"),
                        ("C1=CC=C2C(=C1)C=CN2", "c1ccc2[nH]ccc2c1"),
                        ("CN1C=NC2=C1C(=O)N(C)C(=O)N2C", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
                        ("O=C1C=CC=CN1", "O=c1cccc[nH]1"),
                        ("C1=CC2=CC=CC=CC2=C1", "c1ccc2cccc2cc1")]:
      checkpoint kek
      check canonicalSmiles(kek) == canonicalSmiles(arom)

  test "aromatic systems":
    check aromaticAtoms("C1=CC=CC=C1") == 6
    check aromaticAtoms("C1=CC2=CC=CC=CC2=C1") == 10     # azulene (fused 5+7)
    check aromaticAtoms("[CH+]1C=CC=CC=C1") == 7         # tropylium
    check aromaticAtoms("[CH-]1C=CC=C1") == 5            # cyclopentadienide
    check aromaticAtoms("C1=CC=[O+]C=C1") == 6           # pyrylium
    check aromaticAtoms("[O-][N+]1=CC=CC=C1") == 6       # pyridine N-oxide
    check aromaticAtoms("C[O+]1C=CC=C1") == 5            # O-methylfuranium
    check aromaticAtoms("C[C-]1C=CC=C[N+]1=O") == 6      # exocyclic =O on N+
    check aromaticAtoms("S=C1C=CC=CN1") == 6             # exocyclic C=S
    check aromaticAtoms("C1=CC2=CC=C3C=CC=C4C=CC(=C1)C2=C34") == 16  # pyrene (perimeter rule)
    check aromaticAtoms("C1=CC=C2C=C3C=CC=CC3=CC2=C1") == 14  # anthracene

  test "non-aromatic rings":
    check aromaticAtoms("C1=CC=C1") == 0                 # cyclobutadiene
    check aromaticAtoms("C1=CC=CC=CC=C1") == 0           # cyclooctatetraene
    check aromaticAtoms("O=C1C=CC(=O)C=C1") == 0         # benzoquinone
    check aromaticAtoms("C1CC=CC=C1") == 0               # cyclohexadiene
    check aromaticAtoms("C1=CCC=C1") == 0                # cyclopentadiene
    check aromaticAtoms("B1C=CC=C1") == 0                # borole (4 pi)
    check aromaticAtoms("C=C1C=CC=C1") == 0              # fulvene
    check aromaticAtoms("B1C=CC=CC=C1") == 0             # borepin (boron excluded, as RDKit)
    check aromaticAtoms("O=S1C=CC=C1") == 0              # thiophene S-oxide

  test "biaryl bond stays single":
    let m = parseSmiles("C1=CC=C(C=C1)C1=CC=CC=C1")
    check aromaticAtoms("C1=CC=C(C=C1)C1=CC=CC=C1") == 12
    check m.bonds[m.bondBetween(3, 6)].order == boSingle

  test "fused systems aromatic only as a whole keep a single shared bond":
    # azulene: the 10-electron perimeter is aromatic, the 5- and 7-rings
    # are not, so (as in RDKit) the shared bond is not aromatic
    let m = parseSmiles("c1ccc2cccc2cc1")
    check aromaticAtoms("c1ccc2cccc2cc1") == 10
    var nArom = 0
    for bd in m.bonds:
      if bd.order == boAromatic: inc nArom
    check nArom == 10
    check m.toSmiles == "c1ccc-2cccc2cc1"
    check canonicalSmiles("C1=CC=C2C=CC=C2C=C1") == m.toSmiles

  test "hydrogens are preserved":
    let m = parseSmiles("C1=CC=CN1")
    check m.molecularFormula == "C4H5N"
    check totalHCount(m, 4) == 1
    check parseSmiles("O=C1C=CC=CN1").molecularFormula == "C5H5NO"
