## Converts between SMILES and SD files, for cross-checking against RDKit.
##
##   molconvert sdf2smi   SD file on stdin -> "SMILES<tab>name" lines
##                        ("ERROR<tab>message" if a record can't be read)
##   molconvert sdf2sdf   SD file on stdin -> SD file, rewritten
##   molconvert smi2sdf   SMILES lines on stdin -> SD file (no coordinates)
##   molconvert smi2smi   SMILES lines on stdin -> SMILES lines
## Options: --v3000 writes V3000 molfiles; --addhs makes hydrogens explicit
## atoms before writing. Data items are kept by sdf2sdf; sdf2smi appends
## them as "key=value" fields.

import std/[os, streams, strutils, sequtils]
import chem

let args = commandLineParams()
if args.len == 0 or args[0] notin ["sdf2smi", "sdf2sdf", "smi2sdf", "smi2smi"]:
  quit "usage: molconvert sdf2smi | sdf2sdf | smi2sdf | smi2smi [--v3000] [--addhs]"
let v3000 = "--v3000" in args
let hs = "--addhs" in args
proc prep(m: Mol): Mol = (if hs: m.addHs else: m)
let input = newFileStream(stdin)

case args[0]
of "sdf2smi":
  for r in sdfRecords(input):
    var fields = if r.error.len > 0: @["ERROR", r.error] else: @[r.mol.toSmiles, r.name]
    for k, v in r.props: fields.add k & "=" & v.replace("\n", "\\n")
    echo fields.join("\t")
of "sdf2sdf":
  for r in sdfRecords(input):
    if r.error.len > 0: stderr.writeLine "skipped: ", r.error
    else: stdout.write prep(r.mol).toSdf(r.name, toSeq(r.props.pairs), v3000)
of "smi2smi":
  for line in input.lines:
    let f = line.strip.split(maxsplit = 1)
    if f.len == 0: continue
    try: echo prep(parseSmiles(f[0])).toSmiles, (if f.len > 1: " " & f[1] else: "")
    except ChemError as e: echo "ERROR ", e.msg
else:
  for line in input.lines:
    let f = line.strip.split(maxsplit = 1)
    if f.len == 0: continue
    try:
      stdout.write prep(parseSmiles(f[0])).toSdf(if f.len > 1: f[1] else: "", v3000 = v3000)
    except ChemError as e:
      stderr.writeLine "skipped ", f[0], ": ", e.msg
