"""Cross-check molfile / SD file reading and writing against RDKit.

Usage: python tools/rdkit_molfile_check.py [smiles_file]
(defaults to RDKit's bundled NCI first_5K.smi). Builds examples/molconvert.

Sets: the NCI molecules, plus the stereoisomers RDKit enumerates from those
with potential stereo, each written by RDKit as
  * a V2000 SD file with 2D coordinates (RDKit's layout and wedges),
  * a V3000 SD file with 2D coordinates,
  * a V2000 SD file with 3D coordinates (ETKDG, a sample of the stereo set),
    without and with hydrogen atoms.
For each, our reading must give the same structure as RDKit's reading of the
same file (compared as RDKit canonical isomeric SMILES of our SMILES output),
and the name and data items must survive. Then our rewrite of the file
(V2000 and V3000) must read back in RDKit as the same structure, also with
hydrogens added (`--addhs`: checks their placement keeps the stereo). Finally,
SD files we write from SMILES (no coordinates) must read back in RDKit as
the same constitution, and in our reader as the same structure apart from
double-bond E/Z (tetrahedral stereo travels in the atom parities), and
our SMILES with hydrogens added must read back in RDKit as the input.
"""
import io, os, random, subprocess, sys
from collections import Counter
from rdkit import Chem, RDConfig, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

RDLogger.DisableLog("rdApp.*")
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RDConfig.RDDataDir, "NCI", "first_5K.smi")
exe = os.path.join(root, "examples", "molconvert")
subprocess.run(["nim", "c", "-d:release", "--hints:off", f"--path:{root}/src", "-o:" + exe,
                os.path.join(root, "examples", "molconvert.nim")], check=True)
random.seed(1)


def run(mode, text, *flags):
    return subprocess.run([exe, mode, *flags], input=text, capture_output=True,
                          text=True, check=True).stdout


def canon(smi, isomeric=True):
    m = Chem.MolFromSmiles(smi)
    return None if m is None else Chem.MolToSmiles(m, isomericSmiles=isomeric)


def sdf(mols, v3000=False):
    out = []
    for k, m in enumerate(mols):
        m = Chem.Mol(m)
        m.SetProp("_Name", f"mol{k}")
        m.SetProp("index", str(k))
        m.SetProp("multi line", f"first\nsecond {k}")
        block = Chem.MolToV3KMolBlock(m) if v3000 else Chem.MolToMolBlock(m)
        out.append(block + f"> <index>\n{k}\n\n> <multi line>\nfirst\nsecond {k}\n\n$$$$\n")
    return "".join(out)


def rdkit_read(text):
    return list(Chem.ForwardSDMolSupplier(io.BytesIO(text.encode())))


def check_read(label, mols, v3000=False):
    text = sdf(mols, v3000)
    expected = [Chem.MolToSmiles(m) if m is not None else None for m in rdkit_read(text)]
    rows = run("sdf2smi", text).splitlines()
    stats = Counter(n=len(rows))
    examples = {}
    for k, (row, exp) in enumerate(zip(rows, expected)):
        f = row.split("\t")
        if f[0] == "ERROR":
            stats["read_error"] += 1; examples.setdefault("read_error", (k, f[1])); continue
        if exp is None:
            continue
        if canon(f[0]) != exp:
            stats["structure"] += 1; examples.setdefault("structure", (k, f[0], exp))
        if f[1] != f"mol{k}" or f[2:] != [f"index={k}", f"multi line=first\\nsecond {k}"]:
            stats["props"] += 1; examples.setdefault("props", (k, f[1:]))
    # rewrite with our writer, read back in RDKit
    for flags in [(), ("--v3000",), ("--addhs",)]:
        back = rdkit_read(run("sdf2sdf", text, *flags))
        key = "rewrite" + "".join(flags)
        good = {r.GetProp("_Name"): Chem.MolToSmiles(r) for r in back if r is not None}
        for k, exp in enumerate(expected):
            if exp is not None and good.get(f"mol{k}") != exp:
                stats[key] += 1; examples.setdefault(key, (k, good.get(f"mol{k}"), exp))
    print(f"== {label}: {dict(stats)}")
    for key, ex in examples.items():
        print(f"   {key}: {ex}")


smiles = [l.split()[0] for l in open(path) if l.strip()]
mols = [m for m in (Chem.MolFromSmiles(s) for s in smiles) if m is not None]
opts = StereoEnumerationOptions(onlyUnassigned=True, maxIsomers=4, unique=True)
stereo = []
for m in mols:
    if Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False):
        stereo += [Chem.MolFromSmiles(Chem.MolToSmiles(x)) for x in EnumerateStereoisomers(m, options=opts)]
stereo = [m for m in stereo if m is not None]

check_read("NCI, V2000 2D", mols)
check_read("NCI, V3000 2D", mols, v3000=True)
check_read("stereoisomers, V2000 2D", stereo)
check_read("stereoisomers, V3000 2D", stereo, v3000=True)
sample = random.sample(stereo, 600)
embedded = []
for m in sample:
    mh = Chem.AddHs(m)
    if AllChem.EmbedMolecule(mh, randomSeed=7) == 0:
        embedded.append(mh)
check_read(f"stereoisomers, V2000 3D ({len(embedded)} embedded)", [Chem.RemoveHs(m) for m in embedded])
check_read(f"stereoisomers, V2000 3D with H atoms", embedded)

# writing from SMILES: no coordinates
for label, ms in [("NCI", mols), ("stereoisomers", stereo)]:
    # name each record by its line so skipped SMILES don't shift the rest
    smis = [Chem.MolToSmiles(m) for m in ms]
    text = run("smi2sdf", "\n".join(f"{s} {k}" for k, s in enumerate(smis)))
    ours = {}
    for line in run("sdf2smi", text).splitlines():
        f = line.split("\t")
        if f[0] != "ERROR": ours[int(f[1])] = f[0]
    back = {int(r.GetProp("_Name")): r for r in rdkit_read(text) if r is not None}
    stats = Counter(n=len(smis), not_written=len(smis) - len(ours))
    for k, s in enumerate(smis):
        if k not in ours: continue
        r, o = back.get(k), ours[k]
        if r is None or Chem.MolToSmiles(r, isomericSmiles=False) != canon(s, False):
            stats["rdkit_constitution"] += 1
        # our round trip: same, apart from double-bond E/Z
        t = Chem.MolFromSmiles(s)
        for b in t.GetBonds():
            b.SetStereo(Chem.BondStereo.STEREONONE)
            b.SetBondDir(Chem.BondDir.NONE)
        if canon(o) != canon(Chem.MolToSmiles(t)):
            stats["our_roundtrip"] += 1
            if stats["our_roundtrip"] <= 3: print("   ", s, o)
    print(f"== {label}, written from SMILES (0D): {dict(stats)}")

# hydrogens added to molecules from SMILES
for label, ms in [("NCI", mols), ("stereoisomers", stereo)]:
    smis = [Chem.MolToSmiles(m) for m in ms]
    rows = run("smi2smi", "\n".join(smis), "--addhs").splitlines()
    stats = Counter(n=len(smis))
    for s, row in zip(smis, rows):
        if row.startswith("ERROR"):
            stats["error"] += 1; continue
        if canon(row) != canon(s):
            stats["changed"] += 1
            if stats["changed"] <= 3: print("   ", s, row)
    print(f"== {label}, hydrogens added (SMILES): {dict(stats)}")
