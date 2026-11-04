# chem

Basic cheminformatics for Nim, in the spirit of (a very small part of) RDKit.

```nim
import chem

let m = parseSmiles("O=C(O)c1ccccc1OC(C)=O")
echo m.toSmiles            # CC(=O)Oc1ccccc1C(=O)O   (canonical)
echo m.molecularFormula    # C9H8O4
echo m.molWeight           # 180.159
echo canonicalSmiles("OCC") == canonicalSmiles("C(O)C")   # true
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
  counting). Parsing sanitizes by default, so `C1=CC=CC=C1` and `c1ccccc1`
  give the same molecule and canonical SMILES; impossible aromatic input
  raises `KekulizeError`.
- **Rings**: SSSR (`sssr`), ring bonds/atoms, ring count.
- **Stereochemistry**: tetrahedral centres (`@`/`@@`, including sulfoxides,
  phosphines and bridgehead/aziridine N) and double-bond E/Z, stored as
  cis/trans on the double bond rather than as `/` `\` marks. Parsing removes
  stereo that isn't real (e.g. `CC[C@@](C)(O)CC`), keeping ring cis/trans
  (1,4-disubstituted cyclohexanes) and pseudoasymmetric centres.
- **Canonical ranking** by iterative partition refinement + tie-breaking,
  with stereo labels as invariants, so canonical SMILES is unique for
  stereoisomers too: meso compounds, ring cis/trans and equivalent `/` `\`
  spellings all give one string.
- **Valence model**: implicit hydrogens, charge-aware allowed valences,
  over-valence detection (`ValenceError`, disable with `sanitize = false`).
- **Descriptors**: Hill formula, average MW, monoisotopic mass, atom/heavy
  atom/bond/ring counts, rotatable bonds.
- **Graph API**: build molecules with `addAtom`/`addBond`, `renumberAtoms`
  (stereo-preserving), connected components.

## Modules

| module | contents |
|---|---|
| `chem/mol` | `Mol`, `Atom`, `Bond` types and graph operations |
| `chem/elements` | periodic table: symbols, masses, default valences |
| `chem/valence` | implicit H, valence validation |
| `chem/rings` | SSSR, ring-bond detection, components, ring count |
| `chem/aromaticity` | Kekulization, aromaticity perception |
| `chem/canon` | canonical atom ranks, symmetry classes |
| `chem/stereo` | stereo perception (`assignStereo`) |
| `chem/sanitization` | `sanitize`: kekulize, valence check, aromaticity, stereo |
| `chem/smiles/parser`, `chem/smiles/writer` | SMILES I/O |
| `chem/descriptors` | formula, masses, counts |

## Testing

```sh
nimble test                          # unit tests
python tools/rdkit_crosscheck.py     # compare against RDKit (needs rdkit)
```

On RDKit's bundled 5K NCI set (4,991 molecules RDKit can read, written in
Kekulé form): all formulas, exact masses and aromatic atom counts agree,
every canonical SMILES round-trips to the same structure in RDKit, and
canonical output is the same for randomised aromatic and Kekulé input
SMILES. One molecule (a Mg chelate) is rejected by our valence check.

The script also checks ~4,500 stereoisomers that RDKit enumerates from
the NCI molecules with potential stereo. Canonical output is the same for
every randomised input; stereo counts and round trips match RDKit, except
for 12 cage-aminal isomers (see limitations).

## Known limitations

- **Aromaticity model**: follows RDKit's default model on everything
  tested, but is an independent implementation. Fused ring combinations are
  only considered up to 6 rings at a time, and neutral boron is never
  aromatic (as in RDKit). With `sanitize = false` no kekulization or
  perception happens, and aromatic flags are kept exactly as written.
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
