## Substructure matches, for cross-checking against RDKit.
##
## `smartsmatch patterns.txt < molecules.smi`: patterns are one SMARTS per
## line, molecules SMILES (first field per line). Prints one line per
## molecule: the input SMILES, then for each pattern a tab-separated field
## of matches ("" if none; matches separated by ';', atoms by ','),
## "ERROR" if the SMILES can't be read, or "PATERR" for a pattern that
## doesn't parse.
##
## `smartsmatch --chiral < pairs.txt`: lines of "SMILES<tab>SMARTS"; prints
## 1 or 0 for whether the SMARTS matches with `useChirality`, or "ERROR".

import std/[os, strutils, sequtils]
import chem

if paramCount() != 1:
  quit "usage: smartsmatch patterns.txt < molecules.smi | smartsmatch --chiral < pairs.txt"

if paramStr(1) == "--chiral":
  for line in stdin.lines:
    let f = line.split('\t')
    if f.len != 2:
      echo "ERROR"
      continue
    try:
      echo int(
        parseSmiles(f[0]).hasSubstructMatch(parseSmarts(f[1]), useChirality = true)
      )
    except ChemError:
      echo "ERROR"
  quit()

var patterns: seq[(bool, Query)]
for line in lines(paramStr(1)):
  try:
    patterns.add (true, parseSmarts(line))
  except SmartsError:
    patterns.add (false, Query())

for line in stdin.lines:
  let fields = line.strip.split(maxsplit = 1)
  if fields.len == 0:
    continue
  var row = @[fields[0]]
  try:
    let m = parseSmiles(fields[0])
    for (ok, q) in patterns:
      if not ok:
        row.add "PATERR"
        continue
      row.add m.substructMatches(q, maxMatches = 100000).mapIt(it.join(",")).join(";")
  except ChemError:
    row.add "ERROR"
  echo row.join("\t")
