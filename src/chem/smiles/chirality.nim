## Chirality tags in SMILES written order.
##
## A SMILES `@`/`@@` refers to the neighbours in the order they are written;
## a `Mol` stores it relative to `chiralReference`. Converting either way is
## the same parity flip, plus one RDKit quirk, so the parser and the writer
## both go through `reorderChirality`.

import ../[mol, valence]

func ringDigitChiralityFlip(m: Mol, i: int, ringDigits: int): bool =
  ## RDKit compatibility: for a three-coordinate centre with a lone pair and
  ## only single bonds (amine N, phosphine P, sulfonium S+), RDKit inverts
  ## the meaning of `@`/`@@` once per ring-closure digit written on the
  ## atom (C[N@]1CC1(C)C is read as if the ring partner came last).
  ## Sulfoxides (S=O) and four-coordinate centres follow OpenSMILES.
  ringDigits mod 2 == 1 and m.degree(i) == 3 and totalHCount(m, i) == 0 and
    explicitValence(m, i) == 3

func reorderChirality*(
    m: Mol, i: int, tag: ChiralTag, fromOrder, toOrder: openArray[int], ringDigits: int
): ChiralTag =
  ## `tag` (relative to `fromOrder`) re-expressed relative to `toOrder`,
  ## where one of the two is the SMILES written order of atom `i`, which
  ## carries `ringDigits` ring-closure digits.
  if (permutationParity(fromOrder, toOrder) == 1) xor
      ringDigitChiralityFlip(m, i, ringDigits):
    flip(tag)
  else:
    tag
