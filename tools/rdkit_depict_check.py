"""Cross-check 2D layouts against RDKit.

Usage: python tools/rdkit_depict_check.py [smiles_file] [--sheet FILE.html]
(defaults to RDKit's bundled NCI first_5K.smi). Builds examples/molconvert
and examples/depict.

For the NCI molecules and the stereoisomers RDKit enumerates from them:
  * stereo through the drawing: RDKit reads our SD files (laid out from
    SMILES; wedges and double-bond geometry, as RDKit reads 2D molfiles)
    as the same molecule, stereo included;
  * layout quality, ours and RDKit's (Compute2DCoords), measured the same
    way: spread of bond lengths, molecules with clashing atoms (closer
    than 0.3 bond lengths) and with crossing bonds;
  * timing of the whole SMILES -> SD file run, and of RDKit's layout.
With --sheet, also writes an HTML page with our drawing beside RDKit's for
100 random molecules, to look over.
"""
import os, random, subprocess, sys, time
from collections import Counter
from rdkit import Chem, RDConfig, RDLogger
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

RDLogger.DisableLog("rdApp.*")
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
args = [a for a in sys.argv[1:]]
sheet = None
if "--sheet" in args:
    k = args.index("--sheet")
    sheet = args[k + 1]
    del args[k:k + 2]
path = args[0] if args else os.path.join(RDConfig.RDDataDir, "NCI", "first_5K.smi")
exes = {}
for name in ["molconvert", "depict"]:
    exes[name] = os.path.join(root, "examples", name)
    subprocess.run(["nim", "c", "-d:release", "--hints:off", f"--path:{root}/src", "-o:" + exes[name],
                    os.path.join(root, "examples", name + ".nim")], check=True)
random.seed(1)

smiles = [l.split()[0] for l in open(path) if l.strip()]
mols = [m for m in (Chem.MolFromSmiles(s) for s in smiles) if m is not None]
opts = StereoEnumerationOptions(onlyUnassigned=True, maxIsomers=4, unique=True)
stereo = []
for m in mols:
    if Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False):
        stereo += [Chem.MolFromSmiles(Chem.MolToSmiles(x)) for x in EnumerateStereoisomers(m, options=opts)]


def quality(m):
    """(bond-length spread, clash?, crossings) of a molecule's 2D coordinates."""
    conf = m.GetConformer()
    pos = [conf.GetAtomPosition(i) for i in range(m.GetNumAtoms())]
    lens = [(pos[b.GetBeginAtomIdx()] - pos[b.GetEndAtomIdx()]).Length() for b in m.GetBonds()]
    if not lens:
        return 0.0, False, 0
    mean = sum(lens) / len(lens)
    spread = (sum((x - mean) ** 2 for x in lens) / len(lens)) ** 0.5 / mean
    clash = False
    n = m.GetNumAtoms()
    for i in range(n):
        for j in range(i + 1, n):
            if m.GetBondBetweenAtoms(i, j) is None and (pos[i] - pos[j]).Length() < 0.3 * mean:
                clash = True
                break
        if clash:
            break
    def side(p, a, b):
        return (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x)
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in m.GetBonds()]
    crossings = 0
    for k, (a1, b1) in enumerate(bonds):
        for a2, b2 in bonds[k + 1:]:
            if len({a1, b1, a2, b2}) < 4:
                continue
            p1, p2, p3, p4 = pos[a1], pos[b1], pos[a2], pos[b2]
            if side(p3, p1, p2) * side(p4, p1, p2) < 0 and side(p1, p3, p4) * side(p2, p3, p4) < 0:
                crossings += 1
    return spread, clash, crossings


def summary(qs):
    n = len(qs)
    return (f"bond spread {sum(q[0] for q in qs) / n:.4f}, "
            f"clashes {sum(q[1] for q in qs)} ({100 * sum(q[1] for q in qs) / n:.1f}%), "
            f"with crossings {sum(q[2] > 0 for q in qs)} ({100 * sum(q[2] > 0 for q in qs) / n:.1f}%)")


def canon(m):
    return Chem.MolToSmiles(m)


for label, ms in [("NCI", mols), ("stereoisomers", stereo)]:
    smis = [Chem.MolToSmiles(m) for m in ms]
    lines = [f"{s} {k}" for k, s in enumerate(smis)]
    t = time.perf_counter()
    out = subprocess.run([exes["molconvert"], "smi2sdf"], input="\n".join(lines) + "\n",
                         capture_output=True, text=True, check=True).stdout
    ours_time = time.perf_counter() - t
    stats = Counter(n=len(ms))
    examples = {}
    ours_q = []
    for block in out.split("$$$$\n"):
        if not block.strip():
            continue
        k = int(block.split("\n", 1)[0])
        back = Chem.MolFromMolBlock(block)
        if back is None:
            stats["RDKit can't read"] += 1
            examples.setdefault("RDKit can't read", smis[k])
            continue
        stats["read"] += 1
        if canon(back) != canon(ms[k]):
            flat = Chem.MolToSmiles(back, isomericSmiles=False) == Chem.MolToSmiles(ms[k], isomericSmiles=False)
            key = "stereo differs" if flat else "structure differs"
            stats[key] += 1
            examples.setdefault(key, (smis[k], canon(back)))
        ours_q.append(quality(Chem.MolFromMolBlock(block, sanitize=False, removeHs=False)))
    rd = [Chem.Mol(m) for m in ms]
    t = time.perf_counter()
    for m in rd:
        rdDepictor.Compute2DCoords(m)
    rdkit_time = time.perf_counter() - t
    rdkit_q = [quality(m) for m in rd]
    print(f"== {label}: {dict(stats)}")
    for key, ex in examples.items():
        print(f"   {key}: {ex}")
    print(f"   ours:  {summary(ours_q)}; SMILES -> SD file {ours_time:.2f} s")
    print(f"   RDKit: {summary(rdkit_q)}; Compute2DCoords {rdkit_time:.2f} s")

if sheet:
    sample = random.sample(mols, 100)
    smis = [Chem.MolToSmiles(m) for m in sample]
    ours = subprocess.run([exes["depict"], "lines"], input="\n".join(smis) + "\n",
                          capture_output=True, text=True, check=True).stdout.splitlines()
    rows = []
    for s, m, svg in zip(smis, sample, ours):
        d = rdMolDraw2D.MolDraw2DSVG(300, 220)
        d.DrawMolecule(m)
        d.FinishDrawing()
        rd_svg = d.GetDrawingText().split("?>", 1)[-1]
        ours_svg = svg.split("?>", 1)[-1]
        rows.append(f"<tr><td>{ours_svg}</td><td>{rd_svg}</td><td><code>{s}</code></td></tr>")
    with open(sheet, "w") as f:
        f.write("<!doctype html><meta charset='utf-8'><title>Layouts</title>"
                "<style>td{border:1px solid #ddd;vertical-align:middle} svg{max-width:300px;height:auto}"
                "code{font-size:11px;word-break:break-all}</style>"
                "<table><tr><th>ours</th><th>RDKit</th><th>SMILES</th></tr>" + "".join(rows) + "</table>")
    print("wrote", sheet)
