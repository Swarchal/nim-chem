## SD files: molfiles with named data items, separated by `$$$$`.

import std/[streams, strutils, tables, sequtils]
import ./[mol, molfile]
export tables   # for iterating `props`

type SdfRecord* = object
  mol*: Mol
  name*: string                          ## the molfile's first header line
  props*: OrderedTable[string, string]   ## data items, in file order
  error*: string                         ## why `mol` couldn't be read ("" if it was)

proc parseRecord(lines: seq[string]; sanitize: bool): SdfRecord =
  if lines.len > 0: result.name = lines[0].strip
  var k = 0   # where the data items start
  try:
    (result.mol, k) = readMolLines(lines, 0, sanitize)
  except ChemError as e:
    result.error = e.msg
    while k < lines.len and not lines[k].startsWith(">"): inc k
  # data items: "> <name>" (maybe with other text), value lines, blank line
  while k < lines.len:
    let s = lines[k]
    inc k
    if not s.startsWith(">"): continue
    let open = s.find('<')
    let close = s.find('>', max(open, 1))
    if open < 0 or close < 0: continue
    var value = ""
    while k < lines.len and lines[k].len > 0:
      if value.len > 0: value.add '\n'
      value.add lines[k]
      inc k
    result.props[s[open + 1 ..< close]] = value

template collectRecords(lineSource, sanitize: untyped) =
  ## Yields a record for each `$$$$`-terminated group of lines.
  var lines: seq[string]
  for raw in lineSource:
    var line = raw
    if line.endsWith('\r'): line.setLen line.high
    if line.startsWith("$$$$"):
      yield parseRecord(lines, sanitize)
      lines.setLen 0
    else:
      lines.add line
  if lines.anyIt(it.strip.len > 0):
    yield parseRecord(lines, sanitize)

iterator sdfRecords*(s: Stream; sanitize = true): SdfRecord =
  ## Records of an SD file, read one at a time. A record whose molecule
  ## can't be read is still yielded, with `error` set (and its name and
  ## data items).
  collectRecords(s.lines, sanitize)

iterator sdfRecords*(text: string; sanitize = true): SdfRecord =
  collectRecords(text.splitLines, sanitize)

iterator readSdf*(path: string; sanitize = true): SdfRecord =
  ## Records of the SD file at `path`.
  let s = openFileStream(path)
  defer: s.close()
  for r in sdfRecords(s, sanitize): yield r

proc toSdf*(m: Mol; name = ""; props: openArray[(string, string)] = [];
            v3000 = false): string =
  ## One SD record: the molfile, data items and the `$$$$` line. Raises
  ## `ValueError` for a data item the format can't hold: a name with `<`,
  ## `>` or a line break, or a value with an empty line or a line starting
  ## with `$$$$` (trailing line breaks are dropped).
  result = toMolBlock(m, name, v3000)
  for (key, value) in props:
    if key.contains({'<', '>', '\n', '\r'}):
      raise newException(ValueError, "SD data item name can't be written: " & key)
    let v = value.strip(leading = false, chars = {'\n', '\r'})
    for line in v.splitLines:
      if v.len > 0 and line.strip.len == 0 or line.startsWith("$$$$"):
        raise newException(ValueError, "value of SD data item '" & key &
                           "' has an empty line or a '$$$$' line")
    result &= "> <" & key & ">\n" & v & "\n\n"
  result &= "$$$$\n"

proc toSdf*(r: SdfRecord; v3000 = false): string =
  var props: seq[(string, string)]
  for k, v in r.props: props.add (k, v)
  toSdf(r.mol, r.name, props, v3000)
