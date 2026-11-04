## Shared test helpers.

import std/[random, sequtils]
import chem

proc shuffled*(m: Mol; r: var Rand): Mol =
  ## `m` with its atoms in a random order.
  var order = toSeq(0..<m.numAtoms)
  r.shuffle(order)
  m.renumberAtoms(order)
