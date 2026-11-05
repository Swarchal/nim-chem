## Draws molecules as SVG. Reads SMILES (one per line, first
## whitespace-separated field; the rest of the line is the legend) from
## stdin.
##
##   depict files DIR      one SVG per molecule: DIR/1.svg, DIR/2.svg, ...
##   depict grid FILE      all molecules in one grid SVG (--per-row N, default 4)
##   depict lines          one SVG per input line on stdout, on a single line
##                         ("ERROR<tab>message" for unreadable SMILES)
##
## --highlight SMARTS highlights the first match of SMARTS in each molecule.

import std/[os, strutils]
import chem

let args = commandLineParams()
if args.len == 0 or args[0] notin ["files", "grid", "lines"] or
   (args[0] != "lines" and args.len < 2):
  quit "usage: depict files DIR | grid FILE [--per-row N] | lines  [--highlight SMARTS]"

proc option(name: string): string =
  let k = args.find(name)
  if k >= 0 and k + 1 < args.len: args[k + 1] else: ""

let highlight = if option("--highlight").len > 0: parseSmarts(option("--highlight")) else: Query()
let perRow = if option("--per-row").len > 0: parseInt(option("--per-row")) else: 4

proc match(m: Mol): seq[int] =
  if highlight.numAtoms > 0: m.substructMatch(highlight) else: @[]

proc draw(m: Mol): string =
  m.toSvg(highlightAtoms = m.match, highlightBonds = m.bondsBetween(m.match))

var mols: seq[Mol]
var legends: seq[string]
var matches: seq[seq[int]]
var n = 0
for line in stdin.lines:
  let fields = line.strip.split(maxsplit = 1)
  if fields[0].len == 0: continue
  inc n
  try:
    let m = parseSmiles(fields[0])
    case args[0]
    of "files": writeFile(args[1] / ($n & ".svg"), draw(m))
    of "lines": echo draw(m).replace("\n", "")
    else:
      mols.add m
      matches.add m.match
      legends.add(if fields.len > 1: fields[1] else: fields[0])
  except ChemError as e:
    if args[0] == "lines": echo "ERROR\t", e.msg
    else: stderr.writeLine "line ", n, ": ", e.msg
if args[0] == "grid":
  writeFile(args[1], toSvgGrid(mols, legends, perRow = perRow, highlights = matches))
