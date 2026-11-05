## Reads SMILES (one per line, first whitespace-separated field) from stdin
## and prints the input and its descriptors, tab-separated, in the order of
## the header line (printed first). With --addhs the descriptors are
## computed with hydrogens as atoms (as RDKit does after `AddHs`, which
## changes rotatable bonds, acceptors and Labute ASA).
## Unparseable lines are reported as "ERROR<tab>message".

import std/[os, strutils, strformat, sequtils, math]
import chem

proc descriptors(m: Mol): seq[(string, string)] =
  let crippen = m.crippenContribs
  let r = m.ringCounts
  let q = m.qedProperties
  @[
    ("logp", &"{crippen.mapIt(it.logP).sum:.10f}"),
    ("mr", &"{crippen.mapIt(it.mr).sum:.10f}"),
    ("tpsa", &"{m.tpsa:.10f}"),
    ("labute", &"{m.labuteASA:.10f}"),
    ("hbd", $m.numHDonors),
    ("hba", $m.numHAcceptors),
    ("lipinski_hbd", $m.numLipinskiHDonors),
    ("lipinski_hba", $m.numLipinskiHAcceptors),
    ("heteroatoms", $m.numHeteroatoms),
    ("rotatable", $m.numRotatableBonds),
    ("rotatable_nonstrict", $m.numRotatableBonds(strict = false)),
    ("rings", $r.rings),
    ("aromatic", $r.aromatic),
    ("aliphatic", $r.aliphatic),
    ("saturated", $r.saturated),
    ("heterocycles", $r.heterocycles),
    ("aromatic_carbocycles", $r.aromaticCarbocycles),
    ("aromatic_heterocycles", $r.aromaticHeterocycles),
    ("aliphatic_carbocycles", $r.aliphaticCarbocycles),
    ("aliphatic_heterocycles", $r.aliphaticHeterocycles),
    ("saturated_carbocycles", $r.saturatedCarbocycles),
    ("saturated_heterocycles", $r.saturatedHeterocycles),
    ("fraction_csp3", &"{m.fractionCSP3:.10f}"),
    ("qed", &"{qed(q):.12f}"),
    ("qed_mw", &"{q.mw:.10f}"),
    ("qed_hba", $q.hba),
    ("qed_arom", $q.arom),
    ("qed_alerts", $q.alerts),
  ]

let addH = "--addhs" in commandLineParams()
echo "smiles\t", descriptors(parseSmiles("C")).mapIt(it[0]).join("\t")
for line in stdin.lines:
  let smi = line.strip.split(maxsplit = 1)
  if smi[0].len == 0:
    continue
  try:
    var m = parseSmiles(smi[0])
    if addH:
      m = addHs(m)
    echo smi[0], "\t", descriptors(m).mapIt(it[1]).join("\t")
  except ChemError as e:
    echo &"{smi[0]}\tERROR\t{e.msg}"
