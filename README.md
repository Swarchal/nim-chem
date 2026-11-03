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
- **SMILES writing**, canonical (default) or in input order.
- **Canonical ranking** by iterative partition refinement + tie-breaking.
- **Valence model**: implicit hydrogens, charge-aware allowed valences,
  over-valence detection (`ValenceError`, disable with `sanitize = false`).
- **Descriptors**: Hill formula, average MW, monoisotopic mass, atom/heavy
  atom/bond/ring counts, rotatable bonds.
- **Graph API**: build molecules with `addAtom`/`addBond`, `renumberAtoms`
  (stereo-preserving), ring bonds/atoms, connected components.

## Modules

| module | contents |
|---|---|
| `chem/mol` | `Mol`, `Atom`, `Bond` types and graph operations |
| `chem/elements` | periodic table: symbols, masses, default valences |
| `chem/valence` | implicit H, valence validation |
| `chem/rings` | ring-bond detection, components, ring count |
| `chem/canon` | canonical atom ranks |
| `chem/smiles/parser`, `chem/smiles/writer` | SMILES I/O |
| `chem/descriptors` | formula, masses, counts |

## Testing

```sh
nimble test                          # unit tests
python tools/rdkit_crosscheck.py     # compare against RDKit (needs rdkit)
```

On RDKit's bundled 5K NCI set (4,991 molecules RDKit can read): all formulas
and exact masses agree, every canonical SMILES round-trips to the same
structure in RDKit, and canonical output is the same for randomised input
SMILES. One molecule (a Mg chelate) is rejected by our valence check.

## Known limitations

- **No aromaticity perception or Kekulization**: `C1=CC=CC=C1` and `c1ccccc1`
  are kept as written and give different canonical SMILES. Implicit H on
  aromatic atoms uses a heuristic rather than a Kekulé structure, and
  over-valent aromatic atoms are not detected: `c1ccccc1(C)C` is accepted
  (giving a wrong formula) where RDKit rejects it.
- **Canonical SMILES is not unique for some stereoisomers**: stereo tags are
  carried through correctly but are not used in canonical ranking, and ties
  between symmetric atoms are broken by input order. So the same molecule
  can get more than one canonical string when stereo is involved:
  - meso compounds, e.g. `C[C@H](O)[C@@H](C)O` can come out with the `@`
    and `@@` swapped;
  - cis/trans ring compounds, e.g. `C[C@@H]1CC[C@H](C)CC1`;
  - meaningless tags on non-stereocentres, e.g. `CC[C@@](C)(O)CC`, can come
    out as either `@` or `@@`.

  Don't rely on canonical-string comparison to deduplicate stereoisomers.
- **Double-bond `/` `\` marks are kept as written**, not converted to E/Z,
  so marking different substituents of the same alkene (e.g. `F/C(Cl)=C/F`
  vs `F\C(/Cl)=C/F`) can give different canonical SMILES. On a ring
  closure each mark is read from its own atom, so consistent marks at the
  two ends are opposite symbols (`F/C=C/1.F\1`); matching symbols
  (`F/C=C/1.F/1`) are rejected as contradictory.
- **Canonicalisation is slow for very large molecules**: refinement is
  roughly quadratic on long chains (about 3 s for a 4,000-atom chain).
  Non-canonical output (`toSmiles(canonical = false)`) stays fast.
- Only tetrahedral chirality (`@`, `@@`, `@TH1/2`); `@SP`, `@TB`, `@OH`
  and `@AL` raise a `SmilesError`.
- Average atomic weights are current IUPAC values; RDKit uses some older ones
  (e.g. S 32.067), so MWs can differ by ~0.01 per such atom.
- Exact mass is not corrected for electron mass on charged species (RDKit's
  is; the cross-check script allows for this).
- Valences follow simple default lists, so some valid structures are
  rejected, e.g. Mg with four bonds in chelate complexes.
