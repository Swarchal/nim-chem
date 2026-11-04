## Converts between SMILES and SD files, for cross-checking against RDKit.
##
##   molconvert sdf2smi            SD file on stdin -> "SMILES<tab>name" lines
##                                 ("ERROR<tab>message" if a record can't be read)
##   molconvert sdf2sdf [--v3000]  SD file on stdin -> SD file, rewritten
##   molconvert smi2sdf [--v3000]  SMILES lines on stdin -> SD file (no coordinates)
## Data items are kept by sdf2sdf; sdf2smi appends them as "key=value" fields.

import std/[os, streams, strutils]
import chem

let args = commandLineParams()
if args.len notin 1..2 or args[0] notin ["sdf2smi", "sdf2sdf", "smi2sdf"]:
  quit "usage: molconvert sdf2smi | sdf2sdf [--v3000] | smi2sdf [--v3000]"
let v3000 = "--v3000" in args
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
    else: stdout.write r.toSdf(v3000)
else:
  for line in input.lines:
    let f = line.strip.split(maxsplit = 1)
    if f.len == 0: continue
    try:
      stdout.write parseSmiles(f[0]).toSdf(if f.len > 1: f[1] else: "", v3000 = v3000)
    except ChemError as e:
      stderr.writeLine "skipped ", f[0], ": ", e.msg
