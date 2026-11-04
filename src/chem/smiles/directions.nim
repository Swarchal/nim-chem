## SMILES/SMARTS `/` `\` bond marks and the double-bond stereo they encode,
## shared by the SMILES and SMARTS parsers.

type
  Dir* = enum
    ## a '/' or '\' mark, as read in the direction it was written
    dNone, dUp, dDown

  MarkedStereo* = object
    ## cis/trans read from the marks around one double bond
    bond*: int
    refs*: (int, int)   ## reference neighbours of `bonds[bond].a` and `.b`
    cis*: bool
    conflict*: int      ## a bond whose mark contradicts another one on the
                        ## same end, or -1

func flip(d: Dir): Dir =
  case d
  of dNone: dNone
  of dUp: dDown
  of dDown: dUp

func stereoFromMarks*[G](m: G; isDouble: openArray[bool];
                         marks: openArray[Dir]): seq[MarkedStereo] =
  ## Cis/trans for every double bond of `m` (a `Mol` or `Query`: anything
  ## with `adj` and `bonds[i].a`/`.b`) that has marks at both ends, given
  ## which bonds are double and the mark of each bond as written from
  ## `bonds[i].a`. Each mark is turned into an "outward" direction, read
  ## from the double-bond atom towards its substituent; equal outward marks
  ## at the two ends mean cis (F/C=C\F), different ones trans (F/C=C/F).
  for bi in 0..<m.bonds.len:
    if not isDouble[bi]: continue
    let bd = m.bonds[bi]
    var refs: array[2, (int, Dir)]
    var conflict = -1
    for k, endAtom in [bd.a, bd.b]:
      refs[k] = (-1, dNone)
      for (nb, nbi) in m.adj[endAtom]:
        if nbi == bi: continue
        let d = marks[nbi]
        if d == dNone: continue
        let outward = if m.bonds[nbi].a == endAtom: d else: flip(d)
        if refs[k][0] < 0:
          refs[k] = (nb, outward)
        elif outward == refs[k][1]:
          conflict = nbi
    # Only a conflict if this double bond is marked at both ends: in
    # O=C(/C=C/X)/C=C/Y the carbonyl C's marks belong to the C=C bonds.
    if refs[0][0] >= 0 and refs[1][0] >= 0:
      result.add MarkedStereo(bond: bi, refs: (refs[0][0], refs[1][0]),
                              cis: refs[0][1] == refs[1][1], conflict: conflict)
