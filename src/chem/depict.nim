## 2D coordinates for drawing a molecule (depiction): regular rings,
## 120° zig-zag chains, the right cis/trans on stereo double bonds, and a
## clean-up that flips branches to remove overlaps. Bonds are 1.5 long
## (as RDKit's). The layout depends only on the molecule, not on its atom
## order: it is computed on a canonically numbered copy.
##
## Ring systems (rings sharing atoms) are built ring by ring: the first as
## a regular polygon, each next one by placing its unplaced atoms on an
## arc between placed ones (a regular polygon when it shares one bond;
## spiro rings point away from the system). The ring systems and the
## acyclic atoms joining them form a tree, which is laid out outwards from
## the largest ring system.

import std/[algorithm, math, sequtils]
import ./[mol, rings, canon]

type V = tuple[x, y: float]

const bondLength* = 1.5

func `+`(a, b: V): V = (a.x + b.x, a.y + b.y)
func `-`(a, b: V): V = (a.x - b.x, a.y - b.y)
func `*`(a: V; k: float): V = (a.x * k, a.y * k)
func len(a: V): float = hypot(a.x, a.y)
func cross(a, b: V): float = a.x * b.y - a.y * b.x
func dir(angle: float): V = (cos(angle), sin(angle))
func angleOf(a: V): float = arctan2(a.y, a.x)
func rotate(a: V; angle: float): V =
  let (c, s) = (cos(angle), sin(angle))
  (a.x * c - a.y * s, a.x * s + a.y * c)
func normAngle(a: float): float =
  ## into [0, 2π)
  result = a mod (2 * PI)
  if result < 0: result += 2 * PI

type Layout = object
  m: Mol                  ## the canonically numbered molecule
  pos: seq[V]
  placed: seq[bool]
  queued: seq[bool]       ## already in the expansion queue
  system: seq[int]        ## per atom: its ring system, or -1
  systems: seq[seq[int]]  ## per ring system: its rings (indices into `rings`)
  rings: seq[seq[int]]    ## SSSR, atoms in ring order
  local: seq[seq[V]]      ## per ring system: atom positions in its own frame (by atom)

# ---------------------------------------------------------------------------
# Ring systems, in their own frames

proc placeArc(L: var Layout; pos: var seq[V]; done: var seq[bool];
              a, b: int; run: seq[int]; away: V) =
  ## Places the atoms of `run` (a path from `a` to `b`, both placed) on a
  ## circular arc with bond-length steps, bulging away from `away`.
  let (pa, pb) = (pos[a], pos[b])
  let c = len(pb - pa)
  let segs = float(run.len + 1)
  let mid = (pa + pb) * 0.5
  var nrm: V = (-(pb - pa).y, (pb - pa).x)
  nrm = if len(nrm) > 0: nrm * (1 / len(nrm)) else: (0.0, 1.0)
  if cross(pb - pa, away - pa) > 0: nrm = nrm * -1.0   # bulge away from `away`
  let ratio = c / bondLength
  if c < 1e-9 or ratio >= segs - 1e-9:
    # too far apart (or on top of each other): a straight or bent line
    for k, x in run:
      let t = float(k + 1) / segs
      pos[x] = pa + (pb - pa) * t + (if c < 1e-9: nrm * bondLength * float(k + 1) else: (0.0, 0.0))
      done[x] = true
    return
  # step angle θ: sin(segs θ/2) / sin(θ/2) = c / L, on (0, 2π/segs)
  var (lo, hi) = (1e-9, 2 * PI / segs - 1e-9)
  for _ in 0..80:
    let th = (lo + hi) / 2
    if sin(segs * th / 2) / sin(th / 2) > ratio: lo = th else: hi = th
  let th = (lo + hi) / 2
  let r = bondLength / (2 * sin(th / 2))
  let phi = segs * th / 2
  let centre = mid - nrm * (r * cos(phi))
  # go round from `a` in the direction that passes the bulge side
  let start = angleOf(pa - centre)
  var sign = 1.0
  let probe = centre + dir(start + th) * r
  if (probe - mid).x * nrm.x + (probe - mid).y * nrm.y < (centre + dir(start - th) * r - mid).x * nrm.x +
     (centre + dir(start - th) * r - mid).y * nrm.y:
    sign = -1.0
  for k, x in run:
    pos[x] = centre + dir(start + sign * th * float(k + 1)) * r
    done[x] = true

proc relaxIfClashing(L: Layout; pos: var seq[V]; atoms: seq[int]) =
  ## Cages (cubane, ...) can't be drawn ring by ring without atoms landing
  ## on each other. If that happened, re-lay the system by stress
  ## majorization: distances as close as possible to bond length times the
  ## number of bonds between the atoms, starting from the ring-by-ring
  ## positions.
  let n = atoms.len
  var clash = false
  for p in 0..<n:
    for q in p + 1..<n:
      if L.m.bondBetween(atoms[p], atoms[q]) < 0 and
         len(pos[atoms[p]] - pos[atoms[q]]) < 0.5 * bondLength:
        clash = true
  if not clash: return
  # graph distances within the system
  var local = newSeqWith(L.m.numAtoms, -1)
  for k, x in atoms: local[x] = k
  var dist = newSeqWith(n, newSeqWith(n, high(int)))
  for s in 0..<n:
    dist[s][s] = 0
    var queue = @[s]
    var head = 0
    while head < queue.len:
      let x = queue[head]
      inc head
      for y in L.m.neighbors(atoms[x]):
        let ly = local[y]
        if ly >= 0 and dist[s][ly] == high(int):
          dist[s][ly] = dist[s][x] + 1
          queue.add ly
  var x = atoms.mapIt(pos[it])
  # separate coincident atoms (deterministically)
  for k in 0..<n:
    x[k] = x[k] + dir(2.399963 * float(k)) * (0.05 * bondLength)
  for _ in 1..300:
    for i in 0..<n:
      var acc: V
      var wsum = 0.0
      for j in 0..<n:
        if i == j or dist[i][j] == high(int): continue
        let d = float(dist[i][j]) * bondLength
        let w = 1 / (d * d)
        let delta = x[i] - x[j]
        let l = max(len(delta), 1e-9)
        acc = acc + (x[j] + delta * (d / l)) * w
        wsum += w
      if wsum > 0: x[i] = acc * (1 / wsum)
  for k, a in atoms: pos[a] = x[k]

func ringCisTrans(m: Mol; ring: seq[int]): seq[int] =
  ## Per ring position k: for a stereo double bond from ring[k] to
  ## ring[k+1], +1 if the ring atoms on either side of it are cis, -1 if
  ## trans; else 0.
  let n = ring.len
  result = newSeq[int](n)
  for k in 0..<n:
    let (i, j) = (ring[k], ring[(k + 1) mod n])
    let bi = m.bondBetween(i, j)
    let bd = m.bonds[bi]
    if bd.order != boDouble or bd.stereo == bsNone: continue
    let (a, b) = (ring[(k + n - 1) mod n], ring[(k + 2) mod n])
    let (ri, rj) = if bd.a == i: bd.stereoAtoms else: (bd.stereoAtoms[1], bd.stereoAtoms[0])
    var cis = bd.stereo == bsCis
    if ri != a: cis = not cis
    if rj != b: cis = not cis
    result[k] = if cis: 1 else: -1

proc macrocycle(ring: seq[int]; cisTrans: seq[int]; pos: var seq[V]): bool =
  ## Lays out a large ring with stereo double bonds as a closed walk that
  ## turns at each atom: the same way at both ends of a cis double bond,
  ## opposite ways at a trans one (a regular polygon makes them all cis).
  ## False if the constraints can't be met.
  let n = ring.len
  # turn signs: runs of constrained bonds fix signs up to a flip; flip
  # each run to have as many outward (+1) turns as possible
  var sign = newSeq[int](n)
  var turn = newSeq[float](n)
  let base = 2 * PI / float(n)   # a regular polygon's turn
  let inward = PI / 3            # an inward (crank) turn
  var start = -1
  for k in 0..<n:
    if cisTrans[(k + n - 1) mod n] == 0:
      start = k
      break
  if start < 0: return false   # every bond constrained
  var k = 0
  while k < n:
    let p = (start + k) mod n
    var run = @[p]
    sign[p] = 1
    var q = p
    while cisTrans[q] != 0 and run.len < n:
      let nx = (q + 1) mod n
      sign[nx] = sign[q] * cisTrans[q]
      run.add nx
      q = nx
    if run.mapIt(sign[it]).foldl(a + b) < 0:
      for x in run: sign[x] = -sign[x]
    # the outward turns of a run make up for its inward ones (so that
    # overall the walk follows a circle), and are clear enough to show cis
    if run.len > 1:
      let neg = run.countIt(sign[it] < 0)
      let outward = (float(run.len) * base + float(neg) * inward) / float(run.len - neg)
      for x in run: turn[x] = if sign[x] < 0: -inward else: outward
    k += run.len
  # the atoms outside runs share what is left of the 360°
  var free = 0
  var rest = 2 * PI
  for k in 0..<n:
    if turn[k] == 0: inc free
    else: rest -= turn[k]
  if free == 0: return false
  for k in 0..<n:
    if turn[k] == 0: turn[k] = rest / float(free)
  if turn.anyIt(abs(it) >= PI * 0.95): return false
  var pts = newSeq[V](n)
  var heading = 0.0   # of the bond from atom k-1 to atom k
  for k in 1..<n:
    pts[k] = pts[k - 1] + dir(heading) * bondLength
    heading += turn[k]
  # close the walk: spread the (small) gap along it
  let gap = pts[n - 1] + dir(heading) * bondLength - pts[0]
  for k in 0..<n: pts[k] = pts[k] - gap * (float(k) / float(n))
  # relax: bond lengths, the intended angles (1-3 distances), cis/trans
  # across each double bond (1-4 distance: 2 bond lengths cis, √7 trans),
  # and weakly the rest towards a circle so the ring doesn't crumple
  var target = newSeqWith(n, newSeq[float](n))
  var weight = newSeqWith(n, newSeq[float](n))
  let radius = float(n) * bondLength / (2 * PI)
  for i in 0..<n:
    for j in 0..<n:
      if i == j: continue
      let g = min(abs(i - j), n - abs(i - j))
      target[i][j] = 2 * radius * sin(PI * float(g) / float(n))
      weight[i][j] = 0.05 / (target[i][j] * target[i][j])
  proc setPair(i, j: int; d, w: float) =
    target[i][j] = d
    target[j][i] = d
    weight[i][j] = w / (d * d)
    weight[j][i] = w / (d * d)
  for k in 0..<n:
    let (i, j) = (k, (k + 1) mod n)
    setPair(i, j, bondLength, 1.0)
    setPair(i, (k + 2) mod n, 2 * bondLength * cos(turn[j] / 2), 1.0)
  for k in 0..<n:
    if cisTrans[k] != 0:
      # the 1-4 distance the turns give: same-way turns bring the ends
      # together (cis), opposite ones apart (trans)
      let (ti, tj) = (turn[k], turn[(k + 1) mod n])
      let a: V = dir(PI - ti) * bondLength
      let b: V = (bondLength, 0.0) + dir(tj) * bondLength
      setPair((k + n - 1) mod n, (k + 2) mod n, len(b - a), 2.0)
  proc wrong(pts: seq[V]): int =
    ## double bonds with the wrong cis/trans
    for k in 0..<n:
      if cisTrans[k] != 0:
        let (a, i, j, b) = (pts[(k + n - 1) mod n], pts[k], pts[(k + 1) mod n], pts[(k + 2) mod n])
        if (cross(j - i, a - i) * cross(j - i, b - i) > 0) != (cisTrans[k] > 0): inc result
  let closed = pts
  for _ in 1..300:
    for i in 0..<n:
      var acc: V
      var wsum = 0.0
      for j in 0..<n:
        if i == j: continue
        let delta = pts[i] - pts[j]
        let l = max(len(delta), 1e-9)
        acc = acc + (pts[j] + delta * (target[i][j] / l)) * weight[i][j]
        wsum += weight[i][j]
      pts[i] = acc * (1 / wsum)
  # (relaxing can flip a nearly straight double bond: keep the closed walk then)
  if wrong(pts) > wrong(closed): pts = closed
  for k, x in ring: pos[x] = pts[k]
  true

proc layoutSystem(L: var Layout; s: int) =
  ## Positions for ring system `s` in its own frame.
  let ringIds = L.systems[s]
  let n = L.m.numAtoms
  var pos = newSeq[V](n)
  var done = newSeq[bool](n)
  if ringIds.len == 1 and L.rings[ringIds[0]].len >= 8:
    let ring = L.rings[ringIds[0]]
    let ct = ringCisTrans(L.m, ring)
    if ct.anyIt(it != 0) and macrocycle(ring, ct, pos):
      L.local[s] = pos
      return
  var ringDone = newSeq[bool](ringIds.len)
  # seed: the ring sharing atoms with the most others, then the largest
  proc shared(r1, r2: seq[int]): int = r1.countIt(it in r2)
  var seed = 0
  var best = (-1, -1)
  for k, ri in ringIds:
    let fused = ringIds.countIt(it != ri and shared(L.rings[ri], L.rings[it]) >= 2)
    if (fused, L.rings[ri].len) > best:
      best = (fused, L.rings[ri].len)
      seed = k
  block:
    let ring = L.rings[ringIds[seed]]
    let r = bondLength / (2 * sin(PI / float(ring.len)))
    for k, x in ring:
      pos[x] = dir(PI / 2 + 2 * PI * float(k) / float(ring.len) + PI / float(ring.len)) * r
      done[x] = true
    ringDone[seed] = true
  while not ringDone.allIt(it):
    # next: the ring with the most placed atoms
    var next = -1
    var most = 0
    for k, ri in ringIds:
      if ringDone[k]: continue
      let c = L.rings[ri].countIt(done[it])
      if c > most:
        most = c
        next = k
    if next < 0: break   # (can't happen: the system is connected)
    ringDone[next] = true
    let ring = L.rings[ringIds[next]]
    if ring.allIt(done[it]): continue
    var centroid: V
    var nDone = 0
    for x in 0..<n:
      if done[x]:
        centroid = centroid + pos[x]
        inc nDone
    centroid = centroid * (1 / float(nDone))
    if most == 1:
      # spiro: a regular polygon pointing away from the placed neighbours
      let k0 = ring.find(ring.filterIt(done[it])[0])
      let s0 = ring[k0]
      var away: V
      var cnt = 0
      for nb in L.m.neighbors(s0):
        if done[nb]:
          away = away + (pos[s0] - pos[nb])
          inc cnt
      let d = if len(away) > 1e-9: away * (1 / len(away)) else: (1.0, 0.0)
      let r = bondLength / (2 * sin(PI / float(ring.len)))
      let centre = pos[s0] + d * r
      let start = angleOf(pos[s0] - centre)
      for j in 1..<ring.len:
        let x = ring[(k0 + j) mod ring.len]
        pos[x] = centre + dir(start + 2 * PI * float(j) / float(ring.len)) * r
        done[x] = true
      continue
    # every run of unplaced atoms between two placed ones goes on an arc
    let m = ring.len
    var k = 0
    while k < m and not done[ring[k]]: inc k
    for step in 0..<m:
      let i = (k + step) mod m
      if not done[ring[i]] or done[ring[(i + 1) mod m]]: continue
      var run: seq[int]
      var j = (i + 1) mod m
      while not done[ring[j]]:
        run.add ring[j]
        j = (j + 1) mod m
      placeArc(L, pos, done, ring[i], ring[j], run, centroid)
  relaxIfClashing(L, pos, toSeq(0..<n).filterIt(done[it]))
  L.local[s] = pos

# ---------------------------------------------------------------------------
# Laying out the tree of ring systems and chains

func isLinear(m: Mol; i: int): bool =
  ## sp: a triple bond, or two double bonds
  var doubles = 0
  for (_, bi) in m.adj[i]:
    case m.bonds[bi].order
    of boTriple: return true
    of boDouble: inc doubles
    else: discard
  doubles >= 2

proc placeSystem(L: var Layout; s, atom: int; parent: int) =
  ## Puts ring system `s` in place: `atom` (already positioned) keeps its
  ## position and its outward direction points to `parent` (or, without
  ## one, the system stays in its own frame).
  let local = L.local[s]
  var atoms: seq[int]
  for ri in L.systems[s]:
    for x in L.rings[ri]:
      if x notin atoms: atoms.add x
  var angle = 0.0
  var shift: V = (0.0, 0.0)
  if parent >= 0:
    # outward direction of `atom` in the local frame: the middle of the
    # largest gap between its ring neighbours
    var angs: seq[float]
    for nb in L.m.neighbors(atom):
      if L.system[nb] == s: angs.add normAngle(angleOf(local[nb] - local[atom]))
    angs.sort()
    var bestGap = -1.0
    var outward = 0.0
    for k, a in angs:
      let b = if k + 1 < angs.len: angs[k + 1] else: angs[0] + 2 * PI
      if b - a > bestGap:
        bestGap = b - a
        outward = a + (b - a) / 2
    let want = angleOf(L.pos[parent] - L.pos[atom])
    angle = want - outward
    shift = L.pos[atom] - rotate(local[atom], angle)
  for x in atoms:
    L.pos[x] = rotate(local[x], angle) + shift
    L.placed[x] = true

func side(p, a, b: V): float =
  ## which side of the line a→b `p` is on (sign)
  cross(b - a, p - a)

proc stereoPenalty(L: Layout; i: int; assign: seq[(int, V)]): int =
  ## Stereo double bonds at `i` (with the other end placed) that the
  ## positions in `assign` (for i's unplaced neighbours) get wrong.
  let m = L.m
  proc at(x: int): (bool, V) =
    if L.placed[x]: return (true, L.pos[x])
    for (y, p) in assign:
      if y == x: return (true, p)
    (false, (0.0, 0.0))
  for (j, bi) in m.adj[i]:
    let bd = m.bonds[bi]
    if bd.order != boDouble or bd.stereo == bsNone or not L.placed[j]: continue
    let (ri, rj) = if bd.a == i: bd.stereoAtoms else: (bd.stereoAtoms[1], bd.stereoAtoms[0])
    # any placed substituent on each end, with the cis/trans that implies
    var cis = bd.stereo == bsCis
    var si, sj = -1
    var pi, pj: V
    for (x, _) in m.adj[i]:
      if x == j: continue
      let (ok, p) = at(x)
      if ok:
        si = x
        pi = p
        if x != ri: cis = not cis
        break
    for (y, _) in m.adj[j]:
      if y == i: continue
      let (ok, p) = at(y)
      if ok:
        sj = y
        pj = p
        if y != rj: cis = not cis
        break
    if si < 0 or sj < 0: continue
    let same = side(pi, L.pos[i], L.pos[j]) * side(pj, L.pos[i], L.pos[j]) > 0
    if same != cis: inc result

proc expand(L: var Layout; i, parent: int; queue: var seq[(int, int)]) =
  ## Places the unplaced neighbours of `i` around it.
  let m = L.m
  var unplaced: seq[int]
  var angs: seq[float]
  for nb in m.neighbors(i):
    if L.placed[nb]: angs.add normAngle(angleOf(L.pos[nb] - L.pos[i]))
    elif nb notin unplaced: unplaced.add nb
  if unplaced.len == 0: return
  unplaced.sort()
  let k = angs.len + unplaced.len
  # candidate sets of directions for the unplaced neighbours
  var options: seq[seq[float]]
  if angs.len == 0:
    let step = if unplaced.len == 2 and not isLinear(m, i): 2 * PI / 3
               else: 2 * PI / float(unplaced.len)
    options.add toSeq(0..<unplaced.len).mapIt(-PI / 6 + step * float(it))
  elif angs.len == 1:
    let a = angs[0]
    if k == 2:
      if isLinear(m, i): options.add @[a + PI]
      else: options.add [@[a - 2 * PI / 3], @[a + 2 * PI / 3]]
    elif k == 3:
      options.add @[a - 2 * PI / 3, a + 2 * PI / 3]
      options.add @[a + 2 * PI / 3, a - 2 * PI / 3]
    elif k == 4:
      options.add @[a - PI / 2, a + PI, a + PI / 2]
    else:
      options.add toSeq(1..unplaced.len).mapIt(a + 2 * PI * float(it) / float(k))
  else:
    # spread over the largest free gap
    angs.sort()
    var bestGap = -1.0
    var start = 0.0
    for q, a in angs:
      let b = if q + 1 < angs.len: angs[q + 1] else: angs[0] + 2 * PI
      if b - a > bestGap:
        bestGap = b - a
        start = a
    options.add toSeq(1..unplaced.len).mapIt(start + bestGap * float(it) / float(unplaced.len + 1))
    if unplaced.len == 2: options.add @[options[0][1], options[0][0]]
  # zig-zag preference: a chain continues parallel to the bond before `i`
  var prefer: V = (0.0, 0.0)
  if parent >= 0:
    for gp in m.neighbors(parent):
      if gp != i and L.placed[gp]:
        prefer = L.pos[parent] - L.pos[gp]
        break
  var best = -1
  var bestScore = (high(int), Inf)
  for o, angles in options:
    var assign: seq[(int, V)]
    for q, u in unplaced: assign.add (u, L.pos[i] + dir(angles[q]) * bondLength)
    let pen = stereoPenalty(L, i, assign)
    let d = assign[0][1] - L.pos[i]
    let zig = if len(prefer) > 0: -(d.x * prefer.x + d.y * prefer.y) else: 0.0
    if (pen, zig) < bestScore:
      bestScore = (pen, zig)
      best = o
  for q, u in unplaced:
    L.pos[u] = L.pos[i] + dir(options[best][q]) * bondLength
    L.placed[u] = true
  # ring systems after all of `i`'s neighbours, so their cis/trans can be checked
  for u in unplaced:
    if L.system[u] >= 0:
      placeSystem(L, L.system[u], u, i)
      if stereoPenalty(L, u, @[]) > 0:
        # mirror the ring system in the bond to `i`
        var atoms: seq[int]
        for ri in L.systems[L.system[u]]:
          for x in L.rings[ri]:
            if x notin atoms: atoms.add x
        let old = atoms.mapIt(L.pos[it])
        let (pa, pb) = (L.pos[i], L.pos[u])
        let axis = (pb - pa) * (1 / len(pb - pa))
        for x in atoms:
          let v = L.pos[x] - pa
          L.pos[x] = pa + axis * (2 * (v.x * axis.x + v.y * axis.y)) - v
        if stereoPenalty(L, u, @[]) > 0:
          for k, x in atoms: L.pos[x] = old[k]
      for ri in L.systems[L.system[u]]:
        for x in L.rings[ri]:
          if x != u and not L.queued[x]:
            L.queued[x] = true
            queue.add (x, -1)
    L.queued[u] = true
    queue.add (u, i)

# ---------------------------------------------------------------------------
# Clean-up: flip branches to remove overlaps

func clashScore(m: Mol; pos: seq[V]; atoms: seq[int]; inComp: seq[bool]): float =
  ## How crowded the drawing of `atoms` (a component) is: atoms closer
  ## than 0.9 bond lengths (more the closer), plus crossing bonds.
  const near = 0.9 * bondLength
  for p in 0..<atoms.len:
    for q in p + 1..<atoms.len:
      let (i, j) = (atoms[p], atoms[q])
      let d = len(pos[i] - pos[j])
      if d < near and m.bondBetween(i, j) < 0: result += (bondLength / max(d, 0.05)) ^ 2
  # bond crossings
  for b1 in 0..<m.numBonds:
    let (a1, c1) = (m.bonds[b1].a, m.bonds[b1].b)
    if not inComp[a1]: continue
    for b2 in b1 + 1..<m.numBonds:
      let (a2, c2) = (m.bonds[b2].a, m.bonds[b2].b)
      if a2 == a1 or a2 == c1 or c2 == a1 or c2 == c1: continue
      let (p1, p2, p3, p4) = (pos[a1], pos[c1], pos[a2], pos[c2])
      if side(p3, p1, p2) * side(p4, p1, p2) < 0 and side(p1, p3, p4) * side(p2, p3, p4) < 0:
        result += 1

func movedScore(m: Mol; pos: seq[V]; atoms: seq[int]; inComp, moved: seq[bool]): float =
  ## The part of `clashScore` that moving the `moved` atoms rigidly can
  ## change: pairs of a moved and an unmoved atom, and pairs of bonds not
  ## both inside the moved part.
  const near = 0.9 * bondLength
  for i in atoms:
    if not moved[i]: continue
    for j in atoms:
      if moved[j]: continue
      let d = len(pos[i] - pos[j])
      if d < near and m.bondBetween(i, j) < 0: result += (bondLength / max(d, 0.05)) ^ 2
  for b1 in 0..<m.numBonds:
    let (a1, c1) = (m.bonds[b1].a, m.bonds[b1].b)
    if not inComp[a1] or not (moved[a1] or moved[c1]): continue
    let rigid1 = moved[a1] and moved[c1]
    for b2 in 0..<m.numBonds:
      if b2 == b1: continue
      let (a2, c2) = (m.bonds[b2].a, m.bonds[b2].b)
      if not inComp[a2]: continue
      let touched2 = moved[a2] or moved[c2]
      if (rigid1 and moved[a2] and moved[c2]) or (touched2 and b2 < b1): continue
      if a2 == a1 or a2 == c1 or c2 == a1 or c2 == c1: continue
      let (p1, p2, p3, p4) = (pos[a1], pos[c1], pos[a2], pos[c2])
      if side(p3, p1, p2) * side(p4, p1, p2) < 0 and side(p1, p3, p4) * side(p2, p3, p4) < 0:
        result += 1

proc stereoOk(L: Layout): bool =
  ## Every stereo double bond drawn as stored, with each substituent
  ## clearly off the bond's line (else readers can't tell cis from trans).
  for bi, bd in L.m.bonds:
    if bd.order != boDouble or bd.stereo == bsNone: continue
    let (pa, pb) = (L.pos[bd.a], L.pos[bd.b])
    let l = len(pb - pa)
    for (x, other) in [(bd.a, bd.b), (bd.b, bd.a)]:
      for nb in L.m.neighbors(x):
        if nb != other and abs(side(L.pos[nb], pa, pb)) < 0.2 * bondLength * l: return false
    let (ra, rb) = bd.stereoAtoms
    let same = side(L.pos[ra], pa, pb) * side(L.pos[rb], pa, pb) > 0
    if same != (bd.stereo == bsCis): return false
  true

proc cleanUp(L: var Layout; atoms: seq[int]) =
  let m = L.m
  if atoms.len > 250: return   # (quadratic; large molecules keep the plain layout)
  var inComp = newSeq[bool](m.numAtoms)
  for x in atoms: inComp[x] = true
  var score = clashScore(m, L.pos, atoms, inComp)
  if score == 0: return
  let inRing = ringBonds(m)
  let wasOk = stereoOk(L)
  for _ in 1..3:
    var improved = false
    for bi, bd in m.bonds:
      if inRing[bi] or bd.order != boSingle or not inComp[bd.a]: continue
      # the side of `b` (the smaller one is moved)
      var sideB = newSeq[bool](m.numAtoms)
      var stack = @[bd.b]
      sideB[bd.b] = true
      while stack.len > 0:
        let x = stack.pop
        for y in m.neighbors(x):
          if not sideB[y] and not (x == bd.b and y == bd.a):
            sideB[y] = true
            stack.add y
      var moved = toSeq(0..<m.numAtoms).filterIt(sideB[it])
      var pivot = bd.a
      if moved.len * 2 > atoms.len:
        moved = atoms.filterIt(not sideB[it])
        pivot = bd.b
      var isMoved = newSeq[bool](m.numAtoms)
      for x in moved: isMoved[x] = true
      let before = movedScore(m, L.pos, atoms, inComp, isMoved)
      if before == 0: continue   # nothing this bond can fix
      let (pa, pb) = (L.pos[bd.a], L.pos[bd.b])
      var moves: seq[seq[V]]
      # reflection in the bond's line, and small rotations about the pivot
      block:
        let u = (pb - pa) * (1 / len(pb - pa))
        var refl = L.pos
        for x in moved:
          let v = L.pos[x] - pa
          let along = u * (v.x * u.x + v.y * u.y)
          refl[x] = pa + along * 2 - v
        moves.add refl
      for ang in [PI / 6, -PI / 6, PI / 3, -PI / 3]:
        var rot = L.pos
        for x in moved: rot[x] = L.pos[pivot] + rotate(L.pos[x] - L.pos[pivot], ang)
        moves.add rot
      for k, candidate in moves:
        let s = score - before + movedScore(m, candidate, atoms, inComp, isMoved)
        # rotations bend angles, so they must gain more
        if s < score - (if k == 0: 1e-9 else: 1.0):
          let old = L.pos
          L.pos = candidate
          if wasOk and not stereoOk(L):
            L.pos = old
            continue
          score = s
          improved = true
          break
      if score == 0: return
    if not improved: return

# ---------------------------------------------------------------------------

proc orient(L: var Layout; atoms: seq[int]) =
  ## Centres the component and turns its longest extent horizontal (its
  ## first atom on the left).
  var c: V
  for x in atoms: c = c + L.pos[x]
  c = c * (1 / float(atoms.len))
  var (sxx, sxy, syy) = (0.0, 0.0, 0.0)
  for x in atoms:
    let d = L.pos[x] - c
    sxx += d.x * d.x
    sxy += d.x * d.y
    syy += d.y * d.y
  let angle = 0.5 * arctan2(2 * sxy, sxx - syy)
  for x in atoms: L.pos[x] = rotate(L.pos[x] - c, -angle)
  if atoms.len > 1 and L.pos[atoms[0]].x > 1e-6:
    for x in atoms: L.pos[x] = L.pos[x] * -1.0
  # tidy -0.0 and rounding noise
  for x in atoms:
    if abs(L.pos[x].x) < 1e-12: L.pos[x].x = 0
    if abs(L.pos[x].y) < 1e-12: L.pos[x].y = 0

proc layoutCanonical(m: Mol): seq[V] =
  var L = Layout(m: m, pos: newSeq[V](m.numAtoms), placed: newSeq[bool](m.numAtoms),
                 queued: newSeq[bool](m.numAtoms), system: newSeqWith(m.numAtoms, -1))
  L.rings = sssr(m)
  # ring systems: rings sharing atoms
  var ringSys = newSeqWith(L.rings.len, -1)
  for r in 0..<L.rings.len:
    if ringSys[r] >= 0: continue
    ringSys[r] = L.systems.len
    var members = @[r]
    var k = 0
    while k < members.len:
      let a = members[k]
      inc k
      for b in 0..<L.rings.len:
        if ringSys[b] < 0 and L.rings[a].anyIt(it in L.rings[b]):
          ringSys[b] = L.systems.len
          members.add b
    L.systems.add members
  for s, members in L.systems:
    for r in members:
      for x in L.rings[r]: L.system[x] = s
  L.local = newSeq[seq[V]](L.systems.len)
  for s in 0..<L.systems.len: layoutSystem(L, s)

  let comp = connectedComponents(m)
  var compAtoms = newSeq[seq[int]](componentCount(comp))
  for i in 0..<m.numAtoms: compAtoms[comp[i]].add i
  compAtoms.sort(proc (a, b: seq[int]): int = cmp(a[0], b[0]))
  var right = 0.0
  for atoms in compAtoms:
    # start from the largest ring system, else the first atom
    var start = atoms[0]
    var bestSize = 0
    for x in atoms:
      let s = L.system[x]
      if s >= 0:
        var size = 0
        for r in L.systems[s]: size += L.rings[r].len
        if size > bestSize:
          bestSize = size
          start = x
    var queue: seq[(int, int)]
    L.pos[start] = (0.0, 0.0)
    L.placed[start] = true
    if L.system[start] >= 0:
      placeSystem(L, L.system[start], start, -1)
      for x in atoms:
        if L.placed[x]:
          L.queued[x] = true
          queue.add (x, -1)
    else:
      L.queued[start] = true
      queue.add (start, -1)
    var head = 0
    while head < queue.len:
      let (x, parent) = queue[head]
      inc head
      expand(L, x, parent, queue)
    cleanUp(L, atoms)
    orient(L, atoms)
    # next to the components before
    var (lo, hi) = (Inf, -Inf)
    for x in atoms:
      lo = min(lo, L.pos[x].x)
      hi = max(hi, L.pos[x].x)
    let shift = if atoms == compAtoms[0]: -lo else: right + bondLength * 1.5 - lo
    for x in atoms: L.pos[x].x += shift
    right = hi + shift
  # centre the whole drawing
  if m.numAtoms > 0:
    var (lo, hi) = (Inf, -Inf)
    for p in L.pos:
      lo = min(lo, p.x)
      hi = max(hi, p.x)
    let mid = (lo + hi) / 2
    for p in L.pos.mitems: p.x -= mid
  L.pos

proc canonicalCopy(m: Mol): (Mol, seq[int]) =
  ## `m` numbered by canonical rank, with bonds and neighbour lists in a
  ## canonical order; and each original atom's new index.
  let rank = canonicalRanks(m)
  var order = newSeq[int](m.numAtoms)
  for i, r in rank: order[r] = i
  var c: Mol
  for i in order: c.addAtom m.atoms[i]
  var bonds = toSeq(0..<m.numBonds)
  bonds.sort(proc (x, y: int): int =
    let (a1, b1) = (min(rank[m.bonds[x].a], rank[m.bonds[x].b]), max(rank[m.bonds[x].a], rank[m.bonds[x].b]))
    let (a2, b2) = (min(rank[m.bonds[y].a], rank[m.bonds[y].b]), max(rank[m.bonds[y].a], rank[m.bonds[y].b]))
    cmp((a1, b1), (a2, b2)))
  for bi in bonds:
    let bd = m.bonds[bi]
    # each bond from its lower-ranked atom
    var (a, b) = (rank[bd.a], rank[bd.b])
    var refs = (rank[bd.stereoAtoms[0]], rank[bd.stereoAtoms[1]])
    if a > b:
      swap(a, b)
      refs = (refs[1], refs[0])
    let nb = c.addBond(a, b, bd.order)
    if bd.stereo != bsNone:
      c.bonds[nb].stereo = bd.stereo
      c.bonds[nb].stereoAtoms = refs
  for a in c.adj.mitems: a.sort()
  (c, rank)

proc compute2DCoords*(m: var Mol) =
  ## Sets `m.coords` to a 2D layout (z = 0), replacing any coordinates.
  ## Double-bond cis/trans is drawn as stored; tetrahedral stereo is left
  ## to wedges (see `geometry.chooseWedges`). Stereo double bonds in rings
  ## (8 or more atoms) may come out wrong.
  let (c, rank) = canonicalCopy(m)
  let pos = layoutCanonical(c)
  m.coords = newSeq[Point3](m.numAtoms)
  for i in 0..<m.numAtoms:
    let p = pos[rank[i]]
    m.coords[i] = (p.x, p.y, 0.0)

proc with2DCoords*(m: Mol): Mol =
  ## A copy of `m` with 2D coordinates (see `compute2DCoords`).
  result = m
  compute2DCoords(result)
