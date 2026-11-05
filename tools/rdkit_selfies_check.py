"""Cross-check SELFIES encoding and decoding against the reference `selfies`
package (with RDKit for structures).

Usage: python tools/rdkit_selfies_check.py [smiles_file]
(defaults to RDKit's bundled NCI first_5K.smi). Builds examples/selfiesconv.

For the NCI molecules and the stereoisomers RDKit enumerates from them:
  * translating RDKit's Kekulé SMILES must give exactly `selfies.encoder`'s
    SELFIES;
  * our canonical SELFIES, decoded by `selfies.decoder` and read by RDKit,
    must be the same molecule (stereo included);
  * our canonical SELFIES must be the same for randomised input SMILES;
  * decoding `selfies.encoder`'s SELFIES must give exactly
    `selfies.decoder`'s SMILES.
Then random symbol strings must decode exactly as `selfies.decoder` does
and read as the same molecule as RDKit reads: from the robust alphabet
(where SELFIES guarantees a valid molecule), and with stereo, charged,
isotope, bracket and special symbols added, under each preset of semantic
constraints (where our SMILES reader and RDKit's can disagree on unusual
valences).
"""
import os, random, subprocess, sys
from collections import Counter
import selfies as sf
from rdkit import Chem, RDConfig, RDLogger
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

RDLogger.DisableLog("rdApp.*")
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RDConfig.RDDataDir, "NCI", "first_5K.smi")
exe = os.path.join(root, "examples", "selfiesconv")
subprocess.run(["nim", "c", "-d:release", "--hints:off", f"--path:{root}/src", "-o:" + exe,
                os.path.join(root, "examples", "selfiesconv.nim")], check=True)
random.seed(1)


def run(lines, *args):
    out = subprocess.run([exe, *args], input="\n".join(lines) + "\n",
                         capture_output=True, text=True, check=True).stdout.splitlines()
    assert len(out) == len(lines), (len(out), len(lines))
    return out


_remove_all_h = Chem.RemoveHsParameters()
_remove_all_h.removeDefiningBondStereo = True


def canon(smi):
    """RDKit canonical SMILES, with H atoms RDKit keeps for undefined
    double-bond stereo removed (we don't keep them)."""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return Chem.MolToSmiles(Chem.RemoveHs(m, _remove_all_h))


def try_(f, *a, **k):
    try:
        return f(*a, **k)
    except Exception as e:
        return "ERROR " + type(e).__name__


class Tally:
    """Counts per outcome, with the first example of each."""

    def __init__(self, n):
        self.stats = Counter(n=n)
        self.examples = {}

    def note(self, key, example=None):
        self.stats[key] += 1
        if example is not None:
            self.examples.setdefault(key, example)

    def report(self, label):
        print(f"== {label}: {dict(self.stats)}")
        for key, ex in self.examples.items():
            print(f"   {key}: {ex}")


smiles = [l.split()[0] for l in open(path) if l.strip()]
mols = [m for m in (Chem.MolFromSmiles(s) for s in smiles) if m is not None]
opts = StereoEnumerationOptions(onlyUnassigned=True, maxIsomers=4, unique=True)
stereo = []
for m in mols:
    if Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False):
        # re-read: RDKit's random SMILES of freshly enumerated isomers can
        # disagree with its canonical SMILES on ring stereo
        stereo += [Chem.MolFromSmiles(Chem.MolToSmiles(x)) for x in EnumerateStereoisomers(m, options=opts)]

for label, ms in [("NCI", mols), ("stereoisomers", stereo)]:
    t = Tally(len(ms))
    note = t.note

    kek = []
    for m in ms:
        k = Chem.Mol(m)
        Chem.Kekulize(k, clearAromaticFlags=True)
        kek.append(Chem.MolToSmiles(k, kekuleSmiles=True))
    smis = [Chem.MolToSmiles(m) for m in ms]
    # 1. translation of the same Kekulé SMILES
    ref_enc = [try_(sf.encoder, k) for k in kek]
    for k, ours, ref in zip(kek, run(kek, "translate"), ref_enc):
        if ref.startswith("ERROR") and ours.startswith("ERROR"):
            note("both reject")
        elif ours != ref:
            note("translate", (k, ours, ref))
    # 2. canonical SELFIES, decoded by the reference
    ours_enc = run(smis, "encode")
    for s, enc in zip(smis, ours_enc):
        if enc.startswith("ERROR"):
            note("encode error", (s, enc))
        elif canon(try_(sf.decoder, enc)) != canon(s):
            note("structure", (s, enc, try_(sf.decoder, enc)))
    # 3. canonical: the same for randomised SMILES
    rand = [Chem.MolToSmiles(m, doRandom=True) for m in ms]
    for s, a, b in zip(smis, ours_enc, run(rand, "encode")):
        if a != b and not a.startswith("ERROR"):
            note("not canonical", (s, a, b))
    # 4. decoding the reference's SELFIES
    good = [e for e in ref_enc if not e.startswith("ERROR")]
    for e, ours in zip(good, run(good, "decode")):
        ref = sf.decoder(e)
        if ours != ref:
            note("decode", (e, ours, ref))
    t.report(label)

# random SELFIES
alphabet = sorted(sf.get_semantic_robust_alphabet())
extra = ["[/C]", "[\\C]", "[/N]", "[=N-1]", "[C@@H1]", "[C@H1]", "[C@@]", "[N@]", "[13C]",
         "[2H]", "[H]", "[Na]", "[Fe+2]", "[CH0]", "[NH1]", "[NH3+1]", "[O-1]", "[S+1]",
         "[/\\Ring1]", "[-/Ring1]", "[\\-Ring2]", "[//Ring1]", "[#Ring1]", "[=Ring3]",
         "[Branch3]", "[#Branch3]", "[=Branch3]", "[nop]", "[epsilon]", "[CH4]", "[Cl]",
         "[I]", "[P]", "[=P+1]", "[B-1]", "[Si]", "[Se]", "[#N+1]"]
for preset, flag, pool in [("default", [], "robust"), ("default", [], "all"),
                           ("octet_rule", ["--octet"], "all"),
                           ("hypervalent", ["--hypervalent"], "all")]:
    sf.set_semantic_constraints(preset)
    symbols = alphabet + extra if pool == "all" else sorted(sf.get_semantic_robust_alphabet())
    strings = []
    for _ in range(10000):
        parts = []
        for _ in range(random.choice([1, 1, 1, 2, 3])):
            parts.append("".join(random.choice(symbols) for _ in range(random.randint(1, 60))))
        strings.append(".".join(parts))
    t = Tally(len(strings))
    for s, ours, back in zip(strings, run(strings, "decode", *flag), run(strings, "roundtrip", *flag)):
        ref = try_(sf.decoder, s)
        if ours != ref:
            t.note("decode", (s, ours, ref))
            continue
        # the same molecule: our reading of the SMILES vs RDKit's
        r = canon(ref)
        if r is None and back.startswith("ERROR"):
            t.note("neither reads")
        elif r is None:
            t.note("only we read", (s, ref, back))
        elif back.startswith("ERROR"):
            t.note("only RDKit reads", (s, ref, back, r))
        elif canon(back) != r:
            t.note("structure", (s, ref, back, r))
    t.report(f"random SELFIES, {preset} constraints, {pool} symbols")
sf.set_semantic_constraints("default")
