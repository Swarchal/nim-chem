## Reads SMILES (one per line, first whitespace-separated field) from stdin
## and prints tab-separated: input, canonical SMILES, formula, MW, exact mass,
## number of aromatic atoms, chiral atoms, stereo double bonds and aromatic
## bonds.
## Unparseable lines are reported as "ERROR<tab>message".

import std/[strutils, strformat, sequtils]
import chem

for line in stdin.lines:
  let smi = line.strip.split(maxsplit = 1)
  if smi.len == 0: continue
  try:
    let m = parseSmiles(smi[0])
    let nArom = m.atoms.countIt(it.aromatic)
    let nChiral = m.atoms.countIt(it.chiral != ctNone)
    let nStereoBonds = m.bonds.countIt(it.stereo != bsNone)
    let nAromBonds = m.bonds.countIt(it.order == boAromatic)
    echo &"{smi[0]}\t{m.toSmiles}\t{m.molecularFormula}\t{m.molWeight:.4f}\t{m.exactMass:.5f}\t{nArom}\t{nChiral}\t{nStereoBonds}\t{nAromBonds}"
  except ChemError as e:
    echo &"{smi[0]}\tERROR\t{e.msg}"
