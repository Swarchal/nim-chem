# chem

Basic cheminformatics for Nim, in the spirit of (a very small part of) RDKit.

```nim
import chem

let m = parseSmiles("O=C(O)c1ccccc1OC(C)=O")
echo m.toSmiles            # CC(=O)Oc1ccccc1C(=O)O   (canonical)
echo m.molecularFormula    # C9H8O4
echo m.molWeight           # 180.159
echo canonicalSmiles("OCC") == canonicalSmiles("C(O)C")   # true

let acid = parseSmarts("[CX3](=O)[OX2H1]")
echo m.substructMatches(acid)                 # @[@[1, 0, 2]]  (C, =O, OH)
echo m.hasSubstructMatch(parseSmiles("c1ccccc1OC(C)=O"))   # true

for rec in readSdf("compounds.sdf"):
  if rec.error.len == 0:
    echo rec.name, " ", rec.mol.toSmiles, " ", rec.props.getOrDefault("IC50")
writeFile("aspirin.mol", m.toMolBlock("aspirin"))

let fp = m.morganFingerprint(radius = 2, nBits = 2048)   # ECFP4, RDKit's bits
echo tanimoto(fp, parseSmiles("OC(=O)c1ccccc1O").morganFingerprint)   # 0.448...
```

## Features (v0.1)

- **SMILES parsing** (OpenSMILES subset): organic subset and bracket atoms
  (isotope, `@`/`@@`, H count, charge, atom maps), all bond types incl. `/` `\`,
  branches, ring closures (`%nn`), `.` components. Errors raise `SmilesError`
  with a character position.
- **SMILES writing**, canonical (default) or in input order, aromatic
  (default) or Kekulé (`toSmiles(kekule = true)`).
- **Aromaticity**: Kekulization (`kekulize`, via a blossom matching, so odd
  fused systems like azulene work) and perception (`setAromaticity`, Hückel
  4n+2 over SSSR rings and fused combinations, RDKit-compatible electron
  counting; for fused systems that are only aromatic as a whole, such as
  azulene, only the perimeter bonds become aromatic). Parsing sanitizes by
  default, so `C1=CC=CC=C1` and `c1ccccc1` give the same molecule and
  canonical SMILES; impossible aromatic input raises `KekulizeError`.
- **Rings**: SSSR (`sssr`), symmetrized SSSR (`symmetrizedSssr`: every
  ring that isn't a sum of smaller ones, so cubane has 6), ring
  bonds/atoms, ring count.
- **Substructure search**: SMARTS parsing (`parseSmarts`: atom primitives
  `* a A # D H h R r v X x + - @ @@`, isotopes, element symbols, recursive
  `$(...)`, the `! & , ;` operators, bond primitives `- = # $ : ~ @ / \`)
  and backtracking subgraph matching (`substructMatches`,
  `substructMatch`, `hasSubstructMatch`), with RDKit's semantics: matches
  are uniquified by atom set by default, and `useChirality` also checks
  tetrahedral and cis/trans stereo. A molecule can be used as the query
  (`toQuery`), matching elements, aromaticity, charges and bond orders.
- **Stereochemistry**: tetrahedral centres (`@`/`@@`, including sulfoxides,
  phosphines and bridgehead/aziridine N) and double-bond E/Z, stored as
  cis/trans on the double bond rather than as `/` `\` marks. Parsing removes
  stereo that isn't real (e.g. `CC[C@@](C)(O)CC`), keeping ring cis/trans
  (1,4-disubstituted cyclohexanes) and pseudoasymmetric centres.
- **Canonical ranking** by iterative partition refinement + tie-breaking,
  with stereo labels as invariants, so canonical SMILES is unique for
  stereoisomers too: meso compounds, ring cis/trans and equivalent `/` `\`
  spellings all give one string.
- **Molfiles and SD files**: V2000 and V3000 reading (`parseMolBlock`)
  and writing (`toMolBlock`; V3000 on request or above 999 atoms/bonds),
  with charges, isotopes, atom maps, radicals and explicit valences.
  Stereo is read from 3D coordinates, 2D wedge/hash bonds (and "either"
  marks), or atom parities when there are no coordinates; coordinates are
  kept on the molecule (`coords`) and written back, with wedges chosen for
  each chiral atom in 2D. SD files are read record by record from a string,
  stream or file (`sdfRecords`, `readSdf`), keeping names and data items;
  unreadable records are reported rather than stopping the read.
- **Explicit hydrogens**: the SMILES and molfile readers turn plain
  hydrogen atoms into H counts, as RDKit does (`removeHs = false` keeps
  them), so `[H]C` reads as `C`. Hydrogens that carry information stay:
  isotopes, charged or mapped H, H bonded to H, a dummy or several atoms,
  and an H needed for double-bond cis/trans or for a three-coordinate
  stereocentre. `removeHs` and `addHs` convert either way, keeping stereo;
  with coordinates, `addHs` places the new atoms (in the plane for 2D,
  roughly tetrahedral in 3D).
- **Valence model**: implicit hydrogens, charge-aware allowed valences,
  over-valence detection (`ValenceError`, disable with `sanitize = false`).
- **Descriptors**: Hill formula, average MW, monoisotopic mass (exact
  masses for 3,111 isotope labels), atom/heavy atom/bond/ring counts,
  rotatable bonds.
- **Fingerprints**: Morgan circular fingerprints (`morganFingerprint`,
  folded to a bit vector; `morganCounts`, unfolded counts;
  `morganEnvironments`, which atom and radius gave each ID), ECFP-style or
  FCFP-style (`useFeatures`), with the same IDs and bits as RDKit's Morgan
  generator. Similarity: Tanimoto, Dice, cosine and Tversky on bit vectors,
  Tanimoto and Dice on counts, and `bulkTanimoto`. `BitVect` also converts
  to and from bit strings. Safe to call from threads (`gcsafe`).
- **Graph API**: build molecules with `addAtom`/`addBond`, `renumberAtoms`
  (stereo-preserving), connected components.

## Modules

| module | contents |
|---|---|
| `chem/mol` | `Mol`, `Atom`, `Bond` types and graph operations |
| `chem/elements` | periodic table: symbols, masses, default valences |
| `chem/isotopes` | exact isotope masses (generated by `tools/gen_isotopes.py`) |
| `chem/valence` | implicit H, valence validation |
| `chem/rings` | SSSR, ring-bond detection, components, ring count |
| `chem/aromaticity` | Kekulization, aromaticity perception |
| `chem/canon` | canonical atom ranks, symmetry classes |
| `chem/stereo` | stereo perception (`assignStereo`) |
| `chem/sanitization` | `sanitize`: kekulize, valence check, aromaticity, stereo |
| `chem/query` | `Query` (atom/bond expression graphs), `toQuery` |
| `chem/smarts` | SMARTS parser |
| `chem/substruct` | substructure matching |
| `chem/molfile`, `chem/sdf` | molfile and SD file I/O |
| `chem/geometry` | stereo from coordinates, wedges and atom parities, H placement |
| `chem/hydrogens` | `removeHs`, `addHs` |
| `chem/smiles/parser`, `chem/smiles/writer` | SMILES I/O |
| `chem/descriptors` | formula, masses, counts |
| `chem/fingerprints` | Morgan fingerprints, atom invariants |
| `chem/bitvect` | `BitVect`, similarity measures |

## Testing

```sh
nimble test                          # unit tests
python tools/rdkit_crosscheck.py     # compare against RDKit (needs rdkit)
python tools/rdkit_smarts_check.py   # substructure search against RDKit
python tools/rdkit_molfile_check.py  # molfile / SD file I/O against RDKit
python tools/rdkit_fingerprint_check.py  # Morgan fingerprints against RDKit
```

On RDKit's bundled 5K NCI set (4,991 molecules RDKit can read, written in
Kekulé form): all formulas, exact masses and aromatic atom and bond counts agree,
every canonical SMILES round-trips to the same structure in RDKit, and
canonical output is the same for randomised aromatic and Kekulé input
SMILES. One molecule (a Mg chelate) is rejected by our valence check.

The script also checks ~4,500 stereoisomers that RDKit enumerates from
the NCI molecules with potential stereo. Canonical output is the same for
every randomised input; stereo counts and round trips match RDKit, except
for 12 cage-aminal isomers (see limitations).

The SMARTS check runs the 1,103 SMARTS patterns bundled with RDKit (Crippen
atom types, fragment descriptors, R. Lewis' filters, WEHI PAINS) against
the same 4,991 molecules: all parse, and of the 5.5 million
(molecule, pattern) pairs, 4 give different matches, all from the
molecule rather than the pattern: RDKit charge-separates perchloric acid
when reading it, and treats ferrocene's bonds to Fe differently. Matching
with `useChirality` agrees on 160,808 (stereoisomer, query) pairs apart
from the 12 cage aminals. The whole SMARTS run takes about as long as
RDKit's (roughly 4 s per 1,000 molecules for all 1,103 patterns).

The molfile check has RDKit write the NCI molecules and the stereoisomers
as SD files (V2000 and V3000 with 2D coordinates, and V2000 with 3D
coordinates for a sample of 563 embedded stereoisomers), with names and
multi-line data items. Our reading agrees with RDKit's on every molecule
except the two whose bonds we don't support (a ferrocene drawn with dative
bonds, and the Mg chelate), the 12 cage aminals, 4 skewed cage drawings and
one pseudoasymmetric phosphorus in 3D; names and data items all survive;
and every file we rewrite reads back in RDKit as the same structure as our
reading, with or without hydrogens added. The 3D files are also read with
all their H atoms, with the same results. SD files we write from SMILES
read back in RDKit with the same constitution, and in our reader with the
same stereo apart from double-bond E/Z. SMILES written after `addHs` read
back in RDKit as the input (apart from the cage aminals).

The fingerprint check compares, for the NCI molecules, the stereoisomers
and 3,887 edge cases (every isotope RDKit knows and some it doesn't,
labelled dummy atoms, molecules with over 1,000 atoms matching a feature
pattern), the starting atom invariants (connectivity and feature),
unfolded count fingerprints at radius 0 to 3 (ECFP- and FCFP-style) and
folded 2048- and 1024-bit fingerprints: all are identical to RDKit's, as
are Tanimoto, Dice, cosine and Tversky similarities of 20,000 random pairs.
Fingerprints are faster than RDKit's (ECFP4 for the NCI set: 58 vs 116 ms;
FCFP4: 390 vs 763 ms), and so is similarity: Tanimoto on all 25 million
NCI pairs takes 2.3 s (RDKit 3.4 s), or 1.0 s when built with
`--passC:-mpopcnt` (or `-march=native`) so bit counting uses the CPU
instruction.

## Known limitations

- **Aromaticity model**: follows RDKit's default model on everything
  tested, but is an independent implementation. Fused ring combinations are
  only considered up to 6 rings at a time, and neutral boron is never
  aromatic (as in RDKit). With `sanitize = false` no kekulization or
  perception happens, and aromatic flags are kept exactly as written.
- **Fingerprints**: Morgan only (no atom-pair, topological-torsion, MACCS
  or RDKit path fingerprints), and without chirality (RDKit's
  `useChirality` needs CIP labels).
- **Stereo scope**: tetrahedral centres and double-bond cis/trans only; no
  allenes, atropisomers or `@SP`/`@TB`/`@OH`/`@AL` (these raise
  `SmilesError`). Double bonds in rings smaller than 8 carry no stereo.
- **Stereo choices that differ from RDKit**:
  - cage aminals (hexamethylenetetramine-type) lose their N tags, since the
    cage fixes every N (RDKit keeps them);
  - contradictory `/` `\` marks on one double bond raise `SmilesError`
    (RDKit warns and drops the stereo).
- **RDKit compatibility quirk**: for three-coordinate centres with a lone
  pair and only single bonds (amine N, phosphine P, sulfonium S+), RDKit
  reads `@`/`@@` inverted once per ring-closure digit on the atom. We do
  the same when reading and writing, so SMILES round-trip with RDKit, but
  other toolkits that follow OpenSMILES literally may read those few cases
  differently. Sulfoxides and four-coordinate centres are unaffected.
- **Molfiles without coordinates**: there is no 2D layout yet, so a
  molecule from SMILES is written with zero coordinates. Its tetrahedral
  stereo goes into the atom parities, which this library (and e.g. Open
  Babel) reads back but RDKit ignores; double-bond E/Z can't be written.
  With coordinates (e.g. from a molfile) all stereo is written; if a 2D
  drawing leaves no single bond free to wedge for some chiral atom,
  writing raises `ChemError`. Query features (query atoms such as `A`/`Q`/`L`,
  atom lists, query bond types such as "any") raise `MolfileError`;
  dative bonds (type 9), S-groups and other V3000 collections are not read.
  `toSdf` raises `ValueError` for data items SD files can't hold (a blank
  or `$$$$` line in a value, `<`/`>` in a name).
- **Stereo from 2D drawings** follows RDKit's reading on all but 4 of 3,302
  RDKit-drawn stereoisomers, all bridged cages drawn with overlapping
  bonds, where the plain bonds at a wedged centre lie on one side of it:
  we read the drawing as 3D geometry, RDKit by the order of the bonds
  around the centre. When writing, wedges both readings agree on are
  preferred, so other programs read our drawings as we do.
- **Hydrogens kept where RDKit removes them**: a mapped H (RDKit drops the
  map) and an H on a three-coordinate stereocentre such as `C[S@]([H])=O`
  (RDKit drops the stereo). `addHs` positions are not optimised and can
  clash with nearby atoms.
- **SMARTS scope**: no component-level grouping (`(C).(C)`), reaction
  SMARTS, `@?` or non-tetrahedral chirality. `/` `\` marks only give
  cis/trans to plain `=` query bonds, and `@`/`@@` is checked wherever it
  appears in an atom expression. `R` and `r` use `symmetrizedSssr`, which
  agrees with RDKit's ring counts on the NCI set; in a few cages RDKit's
  ring set is smaller (pentaprismane: RDKit finds only the 4-rings).
  Inside recursive SMARTS, a written H on the first atom (`[$([C@H](C)F)]`)
  is read as it is elsewhere, first in the neighbour order; RDKit reads it
  second there (but first outside `$()`), so such patterns give the
  opposite stereo match to RDKit's.
- **Canonical ties under stereo** are resolved by choosing, at each
  tie-break, the option with the smallest stereo labels; when several
  options tie, each is completed and the smallest full encoding wins. This
  search is capped at 5,000 tie-break trials per molecule, after which the
  first option is taken, so a pathologically symmetric molecule carrying
  many unresolved stereo centres could still depend on atom order.
- **Canonicalisation is slow for very large molecules**: refinement is
  roughly quadratic on long chains (about 3 s for a 4,000-atom chain).
  Non-canonical output (`toSmiles(canonical = false)`) stays fast.
- Average atomic weights are current IUPAC values; RDKit uses some older ones
  (e.g. S 32.067), so MWs can differ by ~0.01 per such atom.
- Exact mass is not corrected for electron mass on charged species (RDKit's
  is; the cross-check script allows for this).
- Valences follow simple default lists, so some valid structures are
  rejected, e.g. Mg with four bonds in chelate complexes.
