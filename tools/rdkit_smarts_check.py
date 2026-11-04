"""Cross-check substructure search against RDKit.

Usage: python tools/rdkit_smarts_check.py [smiles_file]
(defaults to RDKit's bundled NCI first_5K.smi). Builds examples/smartsmatch.

Runs every SMARTS pattern bundled with RDKit (Crippen logP/MR types,
fragment descriptors, R. Lewis' filters, WEHI PAINS) against every molecule
and compares the uniquified matches with RDKit's `GetSubstructMatches`:
the same atom sets, and every reported tuple (query atom -> molecule atom)
a real match.

Then matching with `useChirality` on the stereoisomers RDKit enumerates
from the molecules with potential stereo (as in rdkit_crosscheck.py): each
isomer against every isomer of the same molecule as a whole-molecule query,
and against fragment queries made of each stereocentre with all, or all but
one, of its neighbours. Both toolkits number atoms in SMILES input
order, so the lists are compared directly.
"""
import csv, os, random, subprocess, sys, tempfile
from collections import Counter
from rdkit import Chem, RDConfig, RDLogger
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

RDLogger.DisableLog("rdApp.*")
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data = RDConfig.RDDataDir
path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(data, "NCI", "first_5K.smi")
exe = os.path.join(root, "examples", "smartsmatch")
subprocess.run(["nim", "c", "-d:release", "--hints:off", f"--path:{root}/src", "-o:" + exe,
                os.path.join(root, "examples", "smartsmatch.nim")], check=True)


def patterns():
    out = []
    for line in open(os.path.join(data, "Crippen.txt")):
        if line.strip() and not line.startswith("#"):
            out.append(("crippen", line.split("\t")[1]))
    for line in open(os.path.join(data, "FragmentDescriptors.csv")):
        f = line.rstrip("\n").split("\t")
        if len(f) > 2 and not line.startswith("#"):
            out.append(("fragment", f[2]))
    for line in open(os.path.join(data, "SmartsLib", "RLewis_smarts.txt")):
        if line.strip() and not line.startswith("#"):
            out.append(("lewis", line.split()[0]))
    for row in csv.reader(open(os.path.join(data, "Pains", "wehi_pains.csv"))):
        out.append(("pains", row[0]))
    return [(src, s) for src, s in out if Chem.MolFromSmarts(s) is not None]


pats = patterns()
rd_pats = [Chem.MolFromSmarts(s) for _, s in pats]
smiles = [l.split()[0] for l in open(path) if l.strip()]
smiles = [s for s in smiles if Chem.MolFromSmiles(s) is not None]

with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
    f.write("\n".join(s for _, s in pats) + "\n")
    patfile = f.name
out = subprocess.run([exe, patfile], input="\n".join(smiles), capture_output=True,
                     text=True, check=True).stdout.splitlines()
os.unlink(patfile)

paterr = set()
fails = Counter()
examples = {}
errors = 0
checked = 0
for line in out:
    fields = line.split("\t")
    smi = fields[0]
    if fields[1:] == ["ERROR"]:
        errors += 1
        continue
    m = Chem.MolFromSmiles(smi)
    for k, (field, q) in enumerate(zip(fields[1:], rd_pats)):
        if field == "PATERR":
            paterr.add(k)
            continue
        checked += 1
        ours = sorted(tuple(map(int, x.split(","))) for x in field.split(";")) if field else []
        # uniquified matches keep one tuple per atom set, and which one
        # depends on search order: compare atom sets, and check every tuple
        # is a genuine match
        theirs = m.GetSubstructMatches(q, maxMatches=100000)
        valid = set(m.GetSubstructMatches(q, uniquify=False, maxMatches=100000)) if ours else set()
        if (sorted(sorted(t) for t in ours) != sorted(sorted(t) for t in theirs)
                or any(t not in valid for t in ours)):
            fails[k] += 1
            examples.setdefault(k, (smi, len(ours), len(theirs)))

print(f"{len(pats)} patterns x {len(out)} molecules ({errors} unreadable), {checked} pairs checked")
print(f"patterns we can't parse: {len(paterr)}")
for k in sorted(paterr):
    print(f"  {pats[k][0]}: {pats[k][1]}")
print(f"pairs with different matches: {sum(fails.values())} over {len(fails)} patterns")
for k, n in fails.most_common():
    smi, a, b = examples[k]
    print(f"  {n:5d} {pats[k][0]}: {pats[k][1]}\n        e.g. {smi}: ours {a}, RDKit {b}")


# --- stereo-aware matching
random.seed(1)
opts = StereoEnumerationOptions(onlyUnassigned=True, maxIsomers=4, unique=True)
pairs = []
for s in smiles:
    m = Chem.MolFromSmiles(s)
    if not Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False):
        continue
    # one atom order for all isomers; stereo re-perceived as from SMILES
    isos = [Chem.MolFromSmiles(Chem.MolToSmiles(x, canonical=False))
            for x in EnumerateStereoisomers(m, options=opts)]
    isos = [x for x in isos if x is not None]
    if len(isos) < 2:
        continue
    single = len(Chem.GetMolFrags(m)) == 1
    for t in isos:
        tsmi = Chem.MolToSmiles(t)
        pairs += [(tsmi, Chem.MolToSmiles(q, doRandom=True)) for q in isos]
        for a in t.GetAtoms():
            if a.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
                continue
            nbrs = [n.GetIdx() for n in a.GetNeighbors()]
            for drop in [None] + nbrs:
                use = [a.GetIdx()] + [n for n in nbrs if n != drop]
                for q in isos:
                    pairs.append((tsmi, Chem.MolFragmentToSmiles(
                        q, atomsToUse=use, isomericSmiles=True, canonical=False,
                        rootedAtAtom=random.choice(use) if single else -1)))
out = subprocess.run([exe, "--chiral"], input="\n".join(f"{t}\t{q}" for t, q in pairs),
                     capture_output=True, text=True, check=True).stdout.splitlines()
bad = Counter()
for (t, q), o in zip(pairs, out):
    rd = Chem.MolFromSmiles(t).HasSubstructMatch(Chem.MolFromSmarts(q), useChirality=True)
    if o != str(int(rd)):
        bad[t] += 1
print(f"\nuseChirality: {len(pairs)} (molecule, query) pairs, {sum(bad.values())} differ "
      f"over {len(bad)} molecules")
for t, n in bad.most_common():
    print(f"  {n:5d} {t}")
