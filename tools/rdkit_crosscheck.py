"""Cross-check the nim `chem` library against RDKit.

Usage: python tools/rdkit_crosscheck.py [smiles_file]
(defaults to RDKit's bundled NCI first_5K.smi). Builds examples/smiprops.

Checks per molecule:
  * formula (element counts + charge) and exact mass agree with RDKit
  * our canonical SMILES describes the same structure (RDKit canonical of ours == of input)
  * the number of aromatic atoms matches RDKit's default aromaticity model
  * our canonical SMILES is invariant to RDKit-randomised input SMILES, in
    both aromatic and Kekule form
Average-MW differences are only reported: RDKit uses some older standard
atomic weights (e.g. S 32.067 vs IUPAC 32.06).
"""
import os, re, subprocess, sys
from rdkit import Chem, RDConfig, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcExactMolWt
from rdkit.Chem.Descriptors import MolWt

RDLogger.DisableLog("rdApp.*")
N_RANDOM = 3
ELECTRON_MASS = 0.00054857990946
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RDConfig.RDDataDir, "NCI", "first_5K.smi")
exe = os.path.join(root, "examples", "smiprops")
subprocess.run(["nim", "c", "-d:release", "--hints:off", f"--path:{root}/src", "-o:" + exe,
                os.path.join(root, "examples", "smiprops.nim")], check=True)


def formula_counts(f):
    """Element counts + charge, so Hill-order conventions don't matter."""
    m = re.match(r"^(.*?)([+-]\d*)?$", f)
    counts = {el: int(n or 1) for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", m.group(1))}
    return counts, m.group(2) or ""


def run(smiles):
    """Run smiprops; returns rows of input, canonical, formula, MW, exact mass, #aromatic."""
    out = subprocess.run([exe], input="\n".join(smiles), capture_output=True, text=True, check=True).stdout
    return [l.split("\t") for l in out.splitlines()]


# parse every input with RDKit once
mols = []
for line in open(path):
    f = line.split()
    if f:
        rd = Chem.MolFromSmiles(f[0])
        if rd is not None:
            mols.append((f[0], rd, Chem.MolToSmiles(rd)))

# one batch: original inputs, RDKit canonical, then N_RANDOM random SMILES each
inputs = [s for s, _, _ in mols]
base = [can for _, _, can in mols]
def random_smiles(rd, j):
    """Alternate aromatic and Kekule random SMILES."""
    if j % 2:
        k = Chem.Mol(rd)
        Chem.Kekulize(k, clearAromaticFlags=True)
        return Chem.MolToSmiles(k, doRandom=True, kekuleSmiles=True)
    return Chem.MolToSmiles(rd, doRandom=True)

rand = [random_smiles(rd, j) for _, rd, _ in mols for j in range(N_RANDOM)]
rows = run(inputs + base + rand)
n = len(mols)
in_rows, base_rows, rand_rows = rows[:n], rows[n:2 * n], rows[2 * n:]

stats = dict(n=n, parse_err=0, formula=0, exact_mass=0, aromatic=0, structure=0, invariance=0)
fails = {k: [] for k in stats}
mw_diff = 0
for i, ((smi, rd, rd_can), row) in enumerate(zip(mols, in_rows)):
    if row[1] == "ERROR":
        stats["parse_err"] += 1; fails["parse_err"].append((smi, row[2])); continue
    can, formula, mw, exact, n_arom = row[1], row[2], float(row[3]), float(row[4]), int(row[5])
    rd_arom = sum(a.GetIsAromatic() for a in rd.GetAtoms())
    if n_arom != rd_arom:
        stats["aromatic"] += 1; fails["aromatic"].append((smi, can, rd_can))
    if formula_counts(formula) != formula_counts(CalcMolFormula(rd)):
        stats["formula"] += 1; fails["formula"].append((smi, formula, CalcMolFormula(rd)))
    # RDKit subtracts one electron mass per unit of positive charge; we don't
    if abs(exact - ELECTRON_MASS * Chem.GetFormalCharge(rd) - CalcExactMolWt(rd)) > 1e-4:
        stats["exact_mass"] += 1; fails["exact_mass"].append((smi, exact, CalcExactMolWt(rd)))
    if abs(mw - MolWt(rd)) > 2e-4 * mw:
        mw_diff += 1
    back = Chem.MolFromSmiles(can)
    if back is None or Chem.MolToSmiles(back) != rd_can:
        stats["structure"] += 1; fails["structure"].append((smi, can))
    b = base_rows[i]
    for r in rand_rows[N_RANDOM * i: N_RANDOM * (i + 1)]:
        if b[1] != r[1]:
            stats["invariance"] += 1; fails["invariance"].append((b[0], r[0], b[1], r[1])); break

print(stats)
print(f"(info) average MW differs from RDKit for {mw_diff} molecules (atomic-weight edition)")
for k, v in fails.items():
    for f in v[:5]:
        print(k, *f, sep="\t")
