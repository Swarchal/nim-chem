## chem: basic cheminformatics for Nim.
##
## .. code-block:: nim
##   import chem
##   let mol = parseSmiles("CC(=O)Oc1ccccc1C(=O)O")
##   echo mol.molecularFormula   # C9H8O4
##   echo mol.toSmiles           # canonical SMILES

import chem/[elements, mol, valence, rings, aromaticity, canon, stereo, sanitization, descriptors]
import chem/smiles/[parser, writer]

export elements, mol, valence, rings, aromaticity, canon, stereo, sanitization, descriptors, parser, writer

proc canonicalSmiles*(smiles: string): string =
  ## Convenience: parse a SMILES string and write it back canonically.
  parseSmiles(smiles).toSmiles(canonical = true)
