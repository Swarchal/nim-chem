import std/[unittest, math, sequtils]
import chem

func near(a, b: float): bool = abs(a - b) < 1e-9

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

  test "rotatable bonds, as RDKit (strict by default)":
    check parseSmiles("CCCC").numRotatableBonds == 1
    check parseSmiles("c1ccccc1-c1ccccc1").numRotatableBonds == 1
    check parseSmiles("C1CCCCC1").numRotatableBonds == 0
    check parseSmiles("CC(=O)Oc1ccccc1C(=O)O").numRotatableBonds == 2
    check parseSmiles("CC(=O)Oc1ccccc1C(=O)O").numRotatableBonds(strict = false) == 3
    # tBu and amide C-N bonds only count in the non-strict definition
    check parseSmiles("CC(C)(C)C(=O)NC").numRotatableBonds == 0
    check parseSmiles("CC(C)(C)C(=O)NC").numRotatableBonds(strict = false) == 2
    # as in RDKit, with H atoms methyl bonds stay non-rotatable but OH no
    # longer ends a chain
    check parseSmiles("CCCO").addHs.numRotatableBonds == 2

# reference values from RDKit (rdMolDescriptors, Crippen, QED)
suite "RDKit descriptors":
  test "Crippen logP and MR":
    for (smiles, logP, mr) in [("CC(=O)Oc1ccccc1C(=O)O", 1.3101, 44.71030000000002),
                               ("CC(=O)Nc1ccc(O)cc1", 1.3505999999999998, 42.41050000000002),
                               ("C1CC2CCC1CC2", 2.5866000000000007, 34.68199999999999),
                               ("OCC(F)(F)F", 0.5410000000000001, 13.140800000000002),
                               ("C[C-](C)C", 1.62059, 20.264999999999997),
                               ("c1ccc2[nH]ccc2c1", 2.1678999999999995, 38.29870000000001)]:
      checkpoint smiles
      let m = parseSmiles(smiles)
      check near(m.crippenLogP, logP) and near(m.crippenMR, mr)
      check near(m.crippenContribs.mapIt(it.logP).sum, logP)
      check near(m.addHs.crippenLogP, logP)   # the same with H atoms
    check parseSmiles("CC(=O)Nc1ccc(O)cc1").crippenAtomTypes ==
      @["C1", "C5", "O9", "N4", "C22", "C18", "C18", "C23", "O2", "C18", "C18"]
    # H on N is typed H3, before the general H on heteroatoms (as RDKit's C++)
    check near(parseSmiles("CN").crippenLogP, -0.4251)

  test "TPSA and Labute ASA":
    for (smiles, psa, asa) in [("CC(=O)Oc1ccccc1C(=O)O", 63.60000000000001, 74.75705264447721),
                               ("CC(=O)Nc1ccc(O)cc1", 49.33, 64.66690538413799),
                               ("C1CO1", 12.53, 19.21190721820381),       # O in a 3-ring
                               ("c1ccc2[nH]ccc2c1", 15.79, 53.322228072246475),
                               ("C1CC2CCC1CC2", 0.0, 51.282183521547346)]:
      checkpoint smiles
      let m = parseSmiles(smiles)
      check near(m.tpsa, psa) and near(m.labuteASA, asa)
      check near(m.addHs.tpsa, psa)
    check parseSmiles("*").labuteASA == 0.0
    check parseSmiles("*C").labuteASA == parseSmiles("C").labuteASA
    check parseSmiles("**").labuteASA == 0.0   # (RDKit: NaN, from 0 / 0)

  test "H-bond donors and acceptors":
    let m = parseSmiles("CC(=O)Oc1ccccc1C(=O)O")
    check m.numHDonors == 1 and m.numHAcceptors == 3
    check m.numLipinskiHDonors == 1 and m.numLipinskiHAcceptors == 4
    check parseSmiles("c1ccc2[nH]ccc2c1").numHAcceptors == 0   # 3-connected n
    check parseSmiles("c1ccncc1").numHAcceptors == 1
    check parseSmiles("NCC(=O)N").numHAcceptors == 2           # not the amide N
    check parseSmiles("[NH3+]C").numHDonors == 1
    check parseSmiles("[NH3+]C").addHs.numLipinskiHDonors == 3
    check parseSmiles("CCl.[Na+]").numHeteroatoms == 2
    check parseSmiles("NCCN.NCCN.NCCN").lipinskiViolations == 1   # 12 NH, 6 N
    check parseSmiles("CC(=O)Oc1ccccc1C(=O)O").lipinskiViolations == 0

  test "ring counts":
    let r = parseSmiles("c1ccc2c(c1)CC1CCNC1c1ccoc12").ringCounts
    check r.rings == 4 and r.aromatic == 2 and r.aliphatic == 2
    check r.aromaticCarbocycles == 1 and r.aromaticHeterocycles == 1
    check r.aliphaticCarbocycles == 1 and r.aliphaticHeterocycles == 1
    check r.saturated == 1 and r.saturatedHeterocycles == 1 and r.heterocycles == 2
    # the symmetrized SSSR, as RDKit: 3 rings in bicyclo[2.2.2]octane
    check parseSmiles("C1CC2CCC1CC2").ringCounts.rings == 3
    check parseSmiles("CC").ringCounts == RingCounts()

  test "fraction sp3 carbon":
    check parseSmiles("CC(=O)Oc1ccccc1C(=O)O").fractionCSP3 == 1 / 9
    # carbon with four neighbours, H included (as RDKit)
    check parseSmiles("C[C-](C)C").fractionCSP3 == 0.75
    check parseSmiles("C[CH](C)").fractionCSP3 == 2 / 3
    check parseSmiles("CC#C").fractionCSP3 == 1 / 3
    check parseSmiles("O").fractionCSP3 == 0.0

  test "QED":
    # RDKit's QED doctests (from Bickerton et al.'s own implementation)
    check abs(qed(parseSmiles("N=C(CCSCc1csc(N=C(N)N)n1)NS(N)(=O)=O")) - 0.253) < 1e-3
    check abs(qed(parseSmiles("CNC(=NCCSCc1nc[nH]c1C)NC#N")) - 0.234) < 1e-3
    check abs(qed(parseSmiles("CCCCCNC(=N)NN=Cc1c[nH]c2ccc(CO)cc12")) - 0.234) < 1e-3
    let p = qedProperties(parseSmiles("CC(=O)Oc1ccccc1C(=O)O"))
    check near(p.mw, 180.15899999999996) and near(p.alogp, 1.3101) and near(p.psa, 63.6)
    check (p.hba, p.hbd, p.rotb, p.arom, p.alerts) == (4, 1, 2, 1, 2)
    check near(qed(parseSmiles("CC(=O)Oc1ccccc1C(=O)O")), 0.5501217966938848)
    check near(qed(parseSmiles("C1CC2CCC1CC2")), 0.4493363057408926)
    check near(qed(p), qed(p, qedMeanWeights))
    check qed(p, qedUnitWeights) != qed(p)
