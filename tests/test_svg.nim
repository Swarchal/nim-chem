import std/[unittest, strutils, xmlparser, xmltree, sequtils]
import chem

proc parsed(svg: string): XmlNode = parseXml(svg)

proc count(n: XmlNode; tag: string): int =
  for c in n:
    if c.kind == xnElement:
      if c.tag == tag: inc result
      result += count(c, tag)

proc texts(n: XmlNode): seq[string] =
  ## each <text> element's content, tags removed
  for c in n:
    if c.kind == xnElement:
      if c.tag == "text": result.add c.innerText
      else: result.add texts(c)

suite "SVG drawing":
  test "a valid SVG document with bonds and labels":
    let svg = parseSmiles("CC(=O)Oc1ccccc1C(=O)O").toSvg
    let x = parsed(svg)
    check x.tag == "svg"
    check x.attr("viewBox").len > 0
    # 13 bonds; the 4 double bonds have two lines each
    check count(x, "path") >= 13 + 4
    # (the symbol is centred on the atom; H and charge are separate texts)
    check texts(x).countIt(it == "O") == 4 and texts(x).countIt(it == "H") == 1

  test "labels: hydrogens, charges, isotopes":
    let zwitterion = texts(parsed(parseSmiles("[NH3+]CC(=O)[O-]").toSvg))
    check zwitterion == @["N", "H3", "+", "O", "O", "−"]   # H3N+, on the left
    let iso = texts(parsed(parseSmiles("[13CH3]O").toSvg))
    check "13" in iso and "C" in iso and "H3" in iso
    check texts(parsed(parseSmiles("[Na+].[Cl-]").toSvg)) == @["Na", "+", "Cl", "−"]
    check texts(parsed(parseSmiles("CCC").toSvg)).len == 0   # carbon unlabelled
    check "*" in texts(parsed(parseSmiles("*C").toSvg))

  test "wedges and crossed double bonds":
    # a wedge or hash: a filled path, or short lines
    let chiral = parseSmiles("C[C@H](N)C(=O)O").toSvg
    check "Z\"" in chiral or chiral.count("<path") > parseSmiles("CC(N)C(=O)O").toSvg.count("<path")
    # an unspecified double bond that could have cis/trans is drawn crossed:
    # still two lines, but not parallel to the bond
    check parseSmiles("CC=CC").toSvg.count("<path") == parseSmiles("C/C=C/C").toSvg.count("<path")

  test "highlighting":
    let m = parseSmiles("c1ccccc1O")
    let match = m.substructMatch(parseSmarts("cO"))
    let svg = m.toSvg(highlightAtoms = match, highlightBonds = m.bondsBetween(match))
    check count(parsed(svg), "circle") == 2
    check svg.count("#ff8080") == 3   # 2 atoms, 1 bond

  test "options and fitting":
    let m = parseSmiles("CCCCCCCCCCCCCCCCCCCC")
    let x = parsed(m.toSvg(DrawOptions(width: 200, height: 100)))
    check x.attr("width") == "200.0" and x.attr("height") == "100.0"
    check "<rect" notin m.toSvg(DrawOptions(background: ""))
    check m.toSvg == m.toSvg   # deterministic

  test "grid":
    let mols = ["CCO", "c1ccccc1", "CC(=O)O", "[Na+].[Cl-]", "C1CC1"].mapIt(parseSmiles(it))
    let svg = toSvgGrid(mols, legends = @["ethanol", "benzene", "acetic acid", "salt",
                                          "a very long legend that won't fit in the cell"],
                        perRow = 3)
    let x = parsed(svg)
    check x.attr("width") == "750.0" and x.attr("height") == "400.0"
    check count(x, "g") == 5
    let t = texts(x)
    check "ethanol" in t and "salt" in t
    check t.anyIt(it.endsWith("…"))   # shortened

  test "molecules with coordinates are drawn as they are":
    var m = parseSmiles("CCO")
    m.coords = @[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]
    let svg = m.toSvg
    # a straight line: all y coordinates equal
    let ys = svg.split("M ")[1..^1].mapIt(it.split(",")[1].split(" ")[0])
    check ys.allIt(it == ys[0])

  test "grid with highlights":
    let mols = ["CCO", "CC(=O)O"].mapIt(parseSmiles(it))
    let q = parseSmarts("[OX2H]")
    let svg = toSvgGrid(mols, highlights = mols.mapIt(it.substructMatch(q)))
    check count(parsed(svg), "circle") == 2
