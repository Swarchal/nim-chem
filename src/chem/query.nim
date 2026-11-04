## Query molecules for substructure search.
##
## A `Query` is a graph like a `Mol`, but each atom and bond carries a
## boolean expression of primitives (as in SMARTS) instead of fixed
## properties. Queries come from `smarts.parseSmarts` or from a molecule
## via `toQuery`; `substruct` matches them against molecules.

import ./[mol, valence]

type
  SmartsError* = object of ChemError
    pos*: int ## 0-based character offset where the error was detected

  ExprOp* = enum
    eoPrim
    eoNot
    eoAnd
    eoOr

  Expr*[P] = ref object ## A boolean expression over primitives of type `P`.
    case op*: ExprOp
    of eoPrim: prim*: P
    of eoNot: arg*: Expr[P]
    of eoAnd, eoOr: args*: seq[Expr[P]]

  AtomPrimKind* = enum
    apAny ## `*`
    apAtomicNum ## `#n` and element symbols
    apAromatic ## `a` (and lowercase symbols)
    apAliphatic ## `A` (and uppercase symbols)
    apIsotope ## leading mass number
    apCharge ## `+n` / `-n`
    apTotalH ## `Hn`: implicit + explicit hydrogens, incl. H atoms
    apImplicitH
      ## `hn`: hydrogens not stored as atoms, i.e. `H` minus
      ## H neighbours, as in RDKit (`h` alone: at least one)
    apDegree ## `Dn`: explicit connections
    apConnectivity ## `Xn`: connections including hydrogens
    apValence ## `vn`: total bond order including hydrogens
    apRingCount ## `Rn`: number of SSSR rings containing the atom
    apRingSize ## `rn`: size of the smallest SSSR ring containing it
    apRingBonds ## `xn`: number of ring bonds (`x` alone: at least one)
    apInRing ## `R` / `r` with no number
    apRecursive ## `$(...)`

  AtomPrim* = object
    kind*: AtomPrimKind
    value*: int
    sub*: ref Query ## apRecursive: the query whose first atom is this atom
    subId*: int ## apRecursive: index into the match cache

  BondPrim* = enum
    bpDefault ## unwritten: single or aromatic
    bpSingle
    bpDouble
    bpTriple
    bpQuadruple
    bpAromatic
    bpAny ## `~`
    bpRing ## `@`

  QueryAtom* = object
    expr*: Expr[AtomPrim]
    chiral*: ChiralTag ## `@`/`@@`, relative to `nbrOrder`
    nbrOrder*: seq[int] ## neighbours in written order; -1 = an H count
    ringDigits*: int ## ring-closure digits written on the atom
    atomMap*: int ## `:n` (not used in matching)

  QueryBond* = object
    a*, b*: int
    expr*: Expr[BondPrim]
    stereo*: BondStereo ## from `/` `\` marks (checked with useChirality)
    stereoAtoms*: (int, int) ## reference neighbours of `a` and `b`

  QueryNeeds* = enum
    ## molecule properties a query uses (see `needs`)
    qnHydrogens ## H counts (`H`, `h`, `X`)
    qnRingBonds ## ring-bond flags (`@` bonds, `R`/`r` alone, `x`)
    qnRingSets ## ring counts and sizes (`Rn`, `rn`)
    qnValence ## total valence (`v`)

  Query* = object
    atoms*: seq[QueryAtom]
    bonds*: seq[QueryBond]
    adj*: seq[seq[(int, int)]] ## per atom: (neighbour, bond index)
    numRecursive*: int ## recursive queries in the whole tree

func prim*[P](p: P): Expr[P] =
  Expr[P](op: eoPrim, prim: p)
func notE*[P](e: Expr[P]): Expr[P] =
  Expr[P](op: eoNot, arg: e)
func andE*[P](args: varargs[Expr[P]]): Expr[P] =
  if args.len == 1:
    args[0]
  else:
    Expr[P](op: eoAnd, args: @args)
func orE*[P](args: varargs[Expr[P]]): Expr[P] =
  if args.len == 1:
    args[0]
  else:
    Expr[P](op: eoOr, args: @args)

func atomPrim*(kind: AtomPrimKind, value = 0): Expr[AtomPrim] =
  prim(AtomPrim(kind: kind, value: value))

func needs*(q: Query): set[QueryNeeds] =
  ## The molecule properties `q`'s primitives use, including those of
  ## recursive queries.
  func walk(e: Expr[AtomPrim], acc: var set[QueryNeeds]) =
    case e.op
    of eoPrim:
      case e.prim.kind
      of apTotalH, apImplicitH, apConnectivity:
        acc.incl qnHydrogens
      of apInRing, apRingBonds:
        acc.incl qnRingBonds
      of apRingCount, apRingSize:
        acc.incl qnRingSets
      of apValence:
        acc.incl qnValence
      of apRecursive:
        acc.incl needs(e.prim.sub[])
      else:
        discard
    of eoNot:
      walk(e.arg, acc)
    of eoAnd, eoOr:
      for a in e.args:
        walk(a, acc)
  func walk(e: Expr[BondPrim], acc: var set[QueryNeeds]) =
    case e.op
    of eoPrim:
      if e.prim == bpRing:
        acc.incl qnRingBonds
    of eoNot:
      walk(e.arg, acc)
    of eoAnd, eoOr:
      for a in e.args:
        walk(a, acc)
  for a in q.atoms:
    walk(a.expr, result)
  for b in q.bonds:
    walk(b.expr, result)

func numAtoms*(q: Query): int {.inline.} =
  q.atoms.len
func numBonds*(q: Query): int {.inline.} =
  q.bonds.len

func addAtom*(q: var Query, a: QueryAtom): int {.discardable.} =
  q.atoms.add a
  q.adj.add newSeq[(int, int)]()
  q.atoms.high

func bondBetween*(q: Query, a, b: int): int =
  for (nb, bi) in q.adj[a]:
    if nb == b:
      return bi
  -1

func addBond*(q: var Query, a, b: int, expr: Expr[BondPrim]): int {.discardable.} =
  q.bonds.add QueryBond(a: a, b: b, expr: expr)
  let bi = q.bonds.high
  q.adj[a].add (b, bi)
  q.adj[b].add (a, bi)
  bi

func toQuery*(m: Mol): Query =
  ## A query matching `m` as a substructure the way RDKit matches a
  ## molecule built from SMILES: atoms by element and aromaticity (plus
  ## charge and isotope when they are non-zero), bonds by exact order.
  ## Hydrogen counts are not compared, so `CC=O` matches any aldehyde or
  ## ketone carbonyl. Stereo is kept for `useChirality` matching.
  for i, a in m.atoms:
    var parts = @[
      atomPrim(apAtomicNum, a.element),
      atomPrim(if a.aromatic: apAromatic else: apAliphatic),
    ]
    if a.element == 0:
      parts = @[atomPrim(apAny)]
    if a.charge != 0:
      parts.add atomPrim(apCharge, a.charge)
    if a.isotope != 0:
      parts.add atomPrim(apIsotope, a.isotope)
    result.addAtom QueryAtom(
      expr: andE(parts),
      chiral: a.chiral,
      nbrOrder: chiralReference(m, i),
      atomMap: a.atomMap,
    )
  for bi, bd in m.bonds:
    let p =
      case bd.order
      of boSingle: bpSingle
      of boDouble: bpDouble
      of boTriple: bpTriple
      of boQuadruple: bpQuadruple
      of boAromatic: bpAromatic
      of boNone: bpAny
    result.addBond(bd.a, bd.b, prim(p))
    result.bonds[bi].stereo = bd.stereo
    result.bonds[bi].stereoAtoms = bd.stereoAtoms
