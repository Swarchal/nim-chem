"""Cross-check Morgan fingerprints against RDKit.

Usage: python tools/rdkit_fingerprint_check.py [smiles_file]
(defaults to RDKit's bundled NCI first_5K.smi). Builds examples/fingerprint.

For the NCI molecules and the stereoisomers RDKit enumerates from them,
compares, molecule by molecule:
  * the starting atom invariants (connectivity and feature),
  * unfolded count fingerprints (ECFP and FCFP, radius 0 to 3),
  * folded bit fingerprints (radius 2, 2048 and 1024 bits),
and on edge cases (all isotopes RDKit knows and some it doesn't, labelled
dummy atoms, molecules with over 1,000 atoms matching a feature pattern),
plus similarities between random pairs (Tanimoto, Dice, cosine, Tversky),
computed by RDKit from its fingerprints and by our similarity functions
from ours.
"""
import os, random, subprocess, sys
from collections import Counter
from rdkit import Chem, DataStructs, RDConfig, RDLogger
from rdkit.Chem import rdMolDescriptors as D, rdFingerprintGenerator as G
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

RDLogger.DisableLog("rdApp.*")
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RDConfig.RDDataDir, "NCI", "first_5K.smi")
exe = os.path.join(root, "examples", "fingerprint")
subprocess.run(["nim", "c", "-d:release", "--hints:off", f"--path:{root}/src", "-o:" + exe,
                os.path.join(root, "examples", "fingerprint.nim")], check=True)
random.seed(1)


def run(smis, *args):
    out = subprocess.run([exe, *map(str, args)], input="\n".join(smis) + "\n",
                         capture_output=True, text=True, check=True).stdout.splitlines()
    assert len(out) == len(smis)
    return out


def parse_counts(line):
    return {} if not line else {int(k): int(v) for k, v in (x.split(":") for x in line.split())}


smiles = [l.split()[0] for l in open(path) if l.strip()]
mols = [m for m in (Chem.MolFromSmiles(s) for s in smiles) if m is not None]
opts = StereoEnumerationOptions(onlyUnassigned=True, maxIsomers=4, unique=True)
stereo = []
for m in mols:
    if Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False):
        stereo += [x for x in EnumerateStereoisomers(m, options=opts)]

# edge cases: every isotope RDKit has a mass for, and mass numbers outside
# its table, on an atom bonded to carbon; labelled dummy atoms (attachment
# points); molecules with over 1,000 atoms matching a feature pattern
pt = Chem.GetPeriodicTable()
edge = []
for z in range(1, 119):
    sym = pt.GetElementSymbol(z)
    for a in range(1, 300):
        if pt.GetMassForIsotope(z, a) > 0 or a % 37 == 0:
            edge.append(f"C[{a}{sym}]")
edge += [f"[{a}*]C(=O)O" for a in range(1, 60)] + ["[16*]c1ccccc1", "*CC[2*]"]
edge += [".".join(["c1ccccc1"] * 200), ".".join(["CC(=O)O"] * 400), "C" * 1500]
edge = [m for m in (Chem.MolFromSmiles(x) for x in edge) if m is not None]

for label, ms in [("NCI", mols), ("stereoisomers", stereo), ("edge cases", edge)]:
    smis = [Chem.MolToSmiles(m) for m in ms]
    ms = [Chem.MolFromSmiles(s) for s in smis]   # atoms in the order we see them
    stats = Counter(n=len(smis))
    examples = {}
    def compare(key, ours, theirs, k):
        if ours != theirs:
            stats[key] += 1
            examples.setdefault(key, (smis[k], ours, theirs))
    for features in [False, True]:
        tag = "feature" if features else "connectivity"
        flags = ["--features"] if features else []
        rows = run(smis, "invariants", *flags)
        for k, (m, row) in enumerate(zip(ms, rows)):
            if row.startswith("ERROR"):
                stats["error"] += 1; examples.setdefault("error", (smis[k], row)); continue
            ref = list(D.GetFeatureInvariants(m) if features else D.GetConnectivityInvariants(m))
            compare(f"{tag} invariants", [int(x) for x in row.split()], ref, k)
        for radius in range(4):
            gen = G.GetMorganGenerator(radius=radius, atomInvariantsGenerator=
                                       G.GetMorganFeatureAtomInvGen() if features else None)
            rows = run(smis, "counts", radius, *flags)
            for k, (m, row) in enumerate(zip(ms, rows)):
                if not row.startswith("ERROR"):
                    compare(f"{tag} counts r{radius}", parse_counts(row),
                            gen.GetSparseCountFingerprint(m).GetNonzeroElements(), k)
    for nbits in [2048, 1024]:
        gen = G.GetMorganGenerator(radius=2, fpSize=nbits)
        rows = run(smis, "bits", 2, nbits)
        for k, (m, row) in enumerate(zip(ms, rows)):
            if not row.startswith("ERROR"):
                compare(f"bits {nbits}", [int(x) for x in row.split()],
                        list(gen.GetFingerprint(m).GetOnBits()), k)
    print(f"== {label}: {dict(stats)}")
    for key, ex in examples.items():
        print(f"   {key}: {ex}")

# similarity: random pairs, computed by our similarity functions and RDKit's
smis = [Chem.MolToSmiles(m) for m in mols]
keep = [k for k, r in enumerate(run(smis, "counts", 2)) if not r.startswith("ERROR")]
smis = [smis[k] for k in keep]
mols = [mols[k] for k in keep]
gen = G.GetMorganGenerator(radius=2, fpSize=2048)
bits = [gen.GetFingerprint(m) for m in mols]
counts = [gen.GetSparseCountFingerprint(m) for m in mols]
pairs = [(random.randrange(len(mols)), random.randrange(len(mols))) for _ in range(20000)]
rows = run([f"{smis[i]} {smis[j]}" for i, j in pairs], "similarity")
worst = 0.0
for (i, j), row in zip(pairs, rows):
    ours = [float(x) for x in row.split()]
    theirs = [DataStructs.TanimotoSimilarity(bits[i], bits[j]),
              DataStructs.TanimotoSimilarity(counts[i], counts[j]),
              DataStructs.DiceSimilarity(bits[i], bits[j]),
              DataStructs.DiceSimilarity(counts[i], counts[j]),
              DataStructs.CosineSimilarity(bits[i], bits[j]),
              DataStructs.TverskySimilarity(bits[i], bits[j], 0.3, 0.7)]
    worst = max(worst, *(abs(x - y) for x, y in zip(ours, theirs)))
print(f"== similarity on 20,000 random pairs (Tanimoto and Dice on bits and counts, "
      f"cosine, Tversky): largest difference {worst:.2g}")
