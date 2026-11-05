## SELFIES conversion for lines on stdin, for cross-checking against the
## reference `selfies` package.
##
##   selfiesconv encode     SMILES -> SELFIES (canonical, from `toSelfies`)
##   selfiesconv translate  SMILES -> SELFIES, keeping the SMILES order
##                          (`smilesToSelfies`)
##   selfiesconv decode     SELFIES -> SMILES (`selfiesToSmiles`)
##   selfiesconv roundtrip  SELFIES -> canonical SMILES (`parseSelfies`)
## Options: --octet, --hypervalent (semantic constraints), --loose (no
## strict check when encoding). A line that fails gives "ERROR message".

import std/[os, strutils]
import chem

let args = commandLineParams()
if args.len == 0 or args[0] notin ["encode", "translate", "decode", "roundtrip"]:
  quit "usage: selfiesconv encode|translate|decode|roundtrip [--octet|--hypervalent] [--loose]"
let constraints =
  if "--octet" in args:
    scOctetRule
  elif "--hypervalent" in args:
    scHypervalent
  else:
    scDefault
let strict = "--loose" notin args

for line in stdin.lines:
  let s = line.strip
  try:
    echo (
      case args[0]
      of "encode": parseSmiles(s).toSelfies(strict = strict, constraints = constraints)
      of "translate": smilesToSelfies(s, strict, constraints)
      of "decode": selfiesToSmiles(s, constraints)
      else:
        parseSelfies(s, constraints = constraints).toSmiles
    )
  except ChemError as e:
    echo "ERROR ", e.msg.replace("\n", " ")
