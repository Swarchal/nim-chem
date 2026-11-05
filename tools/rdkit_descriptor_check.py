"""Cross-check descriptors against RDKit.

Usage: python tools/rdkit_descriptor_check.py [smiles_file]
(defaults to RDKit's bundled NCI first_5K.smi). Builds
examples/calcdescriptors.

Every descriptor of examples/calcdescriptors is compared with RDKit's for
each molecule both read: Crippen logP and MR, TPSA, Labute ASA, H-bond
donors and acceptors (RDKit's and Lipinski's), heteroatoms, rotatable
bonds (strict and not), the ring counts, fraction Csp3 and QED with its
properties. Then again with hydrogens as atoms (--addhs / Chem.AddHs).
Reports, per descriptor, how many molecules differ, with examples.
"""
import os, subprocess, sys
from collections import Counter
from rdkit import Chem, RDConfig, RDLogger
from rdkit.Chem import QED, rdMolDescriptors as D

RDLogger.DisableLog("rdApp.*")
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RDConfig.RDDataDir, "NCI", "first_5K.smi")
exe = os.path.join(root, "examples", "calcdescriptors")
subprocess.run(["nim", "c", "-d:release", "--hints:off", f"--path:{root}/src", "-o:" + exe,
                os.path.join(root, "examples", "calcdescriptors.nim")], check=True)

strict, nonstrict = D.NumRotatableBondsOptions.Strict, D.NumRotatableBondsOptions.NonStrict
rdkit = {
    "logp": lambda m: D.CalcCrippenDescriptors(m)[0],
    "mr": lambda m: D.CalcCrippenDescriptors(m)[1],
    "tpsa": D.CalcTPSA,
    "labute": D.CalcLabuteASA,
    "hbd": D.CalcNumHBD,
    "hba": D.CalcNumHBA,
    "lipinski_hbd": D.CalcNumLipinskiHBD,
    "lipinski_hba": D.CalcNumLipinskiHBA,
    "heteroatoms": D.CalcNumHeteroatoms,
    "rotatable": lambda m: D.CalcNumRotatableBonds(m, strict),
    "rotatable_nonstrict": lambda m: D.CalcNumRotatableBonds(m, nonstrict),
    "rings": D.CalcNumRings,
    "aromatic": D.CalcNumAromaticRings,
    "aliphatic": D.CalcNumAliphaticRings,
    "saturated": D.CalcNumSaturatedRings,
    "heterocycles": D.CalcNumHeterocycles,
    "aromatic_carbocycles": D.CalcNumAromaticCarbocycles,
    "aromatic_heterocycles": D.CalcNumAromaticHeterocycles,
    "aliphatic_carbocycles": D.CalcNumAliphaticCarbocycles,
    "aliphatic_heterocycles": D.CalcNumAliphaticHeterocycles,
    "saturated_carbocycles": D.CalcNumSaturatedCarbocycles,
    "saturated_heterocycles": D.CalcNumSaturatedHeterocycles,
    "fraction_csp3": D.CalcFractionCSP3,
    # (`qed` from the properties computed once per molecule, below)
    "qed": lambda m, p: QED.qed(m, qedProperties=p),
    "qed_mw": lambda m, p: p.MW,
    "qed_hba": lambda m, p: p.HBA,
    "qed_arom": lambda m, p: p.AROM,
    "qed_alerts": lambda m, p: p.ALERTS,
}

smiles = [l.split()[0] for l in open(path) if l.strip()]
for label, flags, prep in [("NCI", [], lambda m: m), ("NCI, H as atoms", ["--addhs"], Chem.AddHs)]:
    out = subprocess.run([exe, *flags], input="\n".join(smiles) + "\n", capture_output=True,
                         text=True, check=True).stdout.splitlines()
    header = out[0].split("\t")
    stats, examples = Counter(), {}
    for smi, line in zip(smiles, out[1:]):
        fields = line.split("\t")
        assert fields[0] == smi
        m = Chem.MolFromSmiles(smi)
        if m is None or fields[1] == "ERROR":
            stats["not both read"] += 1
            continue
        stats["compared"] += 1
        m = prep(m)
        props = QED.properties(m)
        for name, value in zip(header[1:], fields[1:]):
            ref = rdkit[name](m, props) if name.startswith("qed") else rdkit[name](m)
            if abs(float(value) - ref) > 1e-6:
                stats[name] += 1
                examples.setdefault(name, (smi, value, ref))
    print(f"== {label}: {dict(stats)}")
    for name, ex in examples.items():
        print(f"   {name}: {ex}")
