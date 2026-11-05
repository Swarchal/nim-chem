## Morgan fingerprints for SMILES on stdin, for cross-checking against RDKit.
##
##   fingerprint counts [radius] [--features]   "id:count id:count ..." per line
##   fingerprint bits [radius] [nBits] [--features]   on-bit indices per line
##   fingerprint invariants [--features]        per-atom starting invariants
##   fingerprint similarity    lines "SMILES SMILES" -> Tanimoto (bits),
##                             Tanimoto (counts), Dice (bits), Dice (counts),
##                             cosine, Tversky 0.3/0.7 (radius 2, 2048 bits)
## A line that can't be read gives "ERROR message".

import std/[os, strutils, tables, algorithm, sequtils]
import chem

let args = commandLineParams()
if args.len == 0 or args[0] notin ["counts", "bits", "invariants", "similarity"]:
  quit "usage: fingerprint counts|bits|invariants|similarity [radius] [nBits] [--features]"
let features = "--features" in args
let nums = args[1 ..^ 1].filterIt(it != "--features").mapIt(parseInt(it))
let radius =
  if nums.len > 0:
    nums[0]
  else:
    2
let nBits =
  if nums.len > 1:
    nums[1]
  else:
    2048

for line in stdin.lines:
  let smi = line.strip.split(maxsplit = 1)
  if smi.len == 0:
    continue
  try:
    let m = parseSmiles(smi[0])
    if args[0] == "similarity":
      let other = parseSmiles(smi[1])
      let (a, b) = (m.morganFingerprint, other.morganFingerprint)
      let (ca, cb) = (m.morganCounts, other.morganCounts)
      echo [
        tanimoto(a, b),
        tanimoto(ca, cb),
        dice(a, b),
        dice(ca, cb),
        cosine(a, b),
        tversky(a, b, 0.3, 0.7),
      ].join(" ")
      continue
    case args[0]
    of "counts":
      var items = toSeq(morganCounts(m, radius, useFeatures = features).pairs)
      items.sort()
      echo items.mapIt($it[0] & ":" & $it[1]).join(" ")
    of "bits":
      echo toSeq(morganFingerprint(m, radius, nBits, useFeatures = features).onBits)
        .join(" ")
    else:
      echo morganInvariants(m, features).join(" ")
  except ChemError as e:
    echo "ERROR ", e.msg
