## Basic molecular descriptors.

import std/[algorithm, sequtils, tables]
import ./[mol, elements, valence, rings]

func hydrogenCount*(m: Mol): int =
  ## Total hydrogens: H atoms in the graph plus implicit/explicit H counts.
  for i, a in m.atoms:
    if a.element == 1: inc result
    result += totalHCount(m, i)

func numHeavyAtoms*(m: Mol): int =
  for a in m.atoms:
    if a.element != 1: inc result

func numAtoms*(m: Mol; includeH: bool): int =
  ## Atom count; with `includeH` implicit and explicit-count H are included.
  result = m.numAtoms
  if includeH:
    for i in 0..<m.numAtoms: result += totalHCount(m, i)

func totalCharge*(m: Mol): int =
  for a in m.atoms: result += a.charge

func molecularFormula*(m: Mol): string =
  ## Molecular formula in Hill order (C, H, then alphabetical; purely
  ## alphabetical if there is no carbon), with a charge suffix such as
  ## "+" or "-2". Isotope labels and wildcard atoms are ignored.
  var counts = initCountTable[string]()
  for a in m.atoms:
    if a.element notin [0, 1]: counts.inc symbol(a.element)
  let h = hydrogenCount(m)
  if h > 0: counts.inc("H", h)
  func part(sym: string): string =
    let n = counts[sym]
    if n == 1: sym else: sym & $n
  let syms = toSeq(counts.keys).sorted
  if "C" in counts:
    result = part("C")
    if "H" in counts: result &= part("H")
    for s in syms:
      if s notin ["C", "H"]: result &= part(s)
  else:
    for s in syms: result &= part(s)
  result &= chargeSuffix(totalCharge(m))

func massSum(m: Mol; elementMass: proc (z: int): float {.noSideEffect, nimcall.}): float =
  ## Sum of atom masses (labelled isotopes use their exact mass) plus attached H.
  for i, a in m.atoms:
    result += (if a.isotope > 0: isotopeMass(a.element, a.isotope)
               else: elementMass(a.element))
    result += float(totalHCount(m, i)) * elementMass(1)

func molWeight*(m: Mol): float =
  ## Average molecular weight (g/mol), including all hydrogens.
  massSum(m, avgMass)

func exactMass*(m: Mol): float =
  ## Monoisotopic mass, using the most abundant isotope of each element
  ## (or the labelled isotope). Electron mass is not corrected for charge.
  massSum(m, monoMass)

func numRotatableBonds*(m: Mol): int =
  ## Simple definition: non-ring single bonds between two atoms that each
  ## have at least one other heavy-atom neighbour.
  let inRing = ringBonds(m)
  proc heavyDegree(i: int): int =
    for nb in m.neighbors(i):
      if m.atoms[nb].element != 1: inc result
  for bi, bd in m.bonds:
    if bd.order != boSingle or inRing[bi]: continue
    if m.atoms[bd.a].element == 1 or m.atoms[bd.b].element == 1: continue
    if heavyDegree(bd.a) > 1 and heavyDegree(bd.b) > 1: inc result
