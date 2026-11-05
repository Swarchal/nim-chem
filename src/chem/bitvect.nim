## Fixed-length bit vectors (folded fingerprints) and similarity measures,
## for bit vectors and for sparse count vectors (unfolded fingerprints,
## `CountTable[uint32]`). Following RDKit, every similarity of two empty
## fingerprints is 0.

import std/[bitops, tables, math]

type BitVect* = object
  size: int
  ones: int             ## number of set bits, kept up to date
  words: seq[uint64]

when defined(c) or defined(cpp):
  # the popcount instruction if the build allows it (`--passC:-mpopcnt`,
  # `-march=native`), else a branch-free fallback, which is faster than
  # the library call `countSetBits` makes then
  {.emit: """/*INCLUDESECTION*/
static inline int chemPopcount64(unsigned long long x) {
#if defined(__POPCNT__) || defined(__aarch64__)
  return __builtin_popcountll(x);
#else
  x = x - ((x >> 1) & 0x5555555555555555ULL);
  x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
  x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0fULL;
  return (int)((x * 0x0101010101010101ULL) >> 56);
#endif
}""".}
  func popcount(x: uint64): int {.importc: "chemPopcount64", nodecl.}
else:
  func popcount(x: uint64): int {.inline.} = countSetBits(x)

func initBitVect*(size: int): BitVect =
  ## An all-zero vector of `size` bits.
  if size < 0: raise newException(ValueError, "negative bit vector size")
  BitVect(size: size, words: newSeq[uint64]((size + 63) div 64))

func len*(v: BitVect): int {.inline.} = v.size

func count*(v: BitVect): int {.inline.} =
  ## Number of set bits.
  v.ones

func checkIndex(v: BitVect; i: int) {.inline.} =
  if i notin 0..<v.size:
    raise newException(IndexDefect, "bit " & $i & " out of range 0.." & $(v.size - 1))

func `[]`*(v: BitVect; i: int): bool =
  v.checkIndex(i)
  v.words[i shr 6].testBit(i and 63)

func `[]=`*(v: var BitVect; i: int; on: bool) =
  if v[i] == on: return
  if on:
    v.words[i shr 6].setBit(i and 63)
    inc v.ones
  else:
    v.words[i shr 6].clearBit(i and 63)
    dec v.ones

func incl*(v: var BitVect; i: int) {.inline.} = v[i] = true

iterator onBits*(v: BitVect): int =
  ## Indices of the set bits, in increasing order.
  for k, w in v.words:
    var x = w
    while x != 0:
      yield k * 64 + countTrailingZeroBits(x)
      x = x and (x - 1)

func checkSizes(a, b: BitVect) =
  if a.size != b.size:
    raise newException(ValueError, "bit vectors differ in length: " & $a.size & " and " & $b.size)

template wordwise(name, op: untyped) =
  func name*(a, b: BitVect): BitVect =
    checkSizes(a, b)
    result = initBitVect(a.size)
    for k in 0..<a.words.len:
      result.words[k] = op(a.words[k], b.words[k])
      result.ones += popcount(result.words[k])

wordwise(`and`, `and`)
wordwise(`or`, `or`)

func commonCount(a, b: BitVect): int =
  ## |a ∧ b|: the one pass over the words each similarity needs.
  checkSizes(a, b)
  {.push boundChecks: off.}
  for k in 0..<a.words.len: result += popcount(a.words[k] and b.words[k])
  {.pop.}

func `$`*(v: BitVect): string =
  ## The bits as '0'/'1' characters, bit 0 first (as RDKit's `ToBitString`).
  result = newString(v.size)
  for i in 0..<v.size: result[i] = if v[i]: '1' else: '0'

func parseBitVect*(s: string): BitVect =
  ## The inverse of `$`.
  result = initBitVect(s.len)
  for i, c in s:
    case c
    of '1': result.incl i
    of '0': discard
    else: raise newException(ValueError, "not a bit string: " & s)

func ratio(num, den: float): float {.inline.} =
  if den == 0: 0.0 else: num / den

func tanimoto*(a, b: BitVect): float =
  ## |a ∧ b| / |a ∨ b| (Jaccard).
  let c = commonCount(a, b)
  ratio(c.float, (a.count + b.count - c).float)

func dice*(a, b: BitVect): float =
  ## 2|a ∧ b| / (|a| + |b|).
  ratio(2.0 * commonCount(a, b).float, (a.count + b.count).float)

func cosine*(a, b: BitVect): float =
  ## |a ∧ b| / sqrt(|a| |b|).
  ratio(commonCount(a, b).float, sqrt((a.count * b.count).float))

func tversky*(a, b: BitVect; alpha, beta: float): float =
  ## |a ∧ b| / (alpha |a \ b| + beta |b \ a| + |a ∧ b|): with alpha =
  ## beta = 1 Tanimoto, with 0.5 Dice; alpha = 1, beta = 0 measures how
  ## much of `a` is in `b`.
  let c = commonCount(a, b).float
  ratio(c, alpha * (a.count.float - c) + beta * (b.count.float - c) + c)

func countSums(a, b: CountTable[uint32]): (int, int, int) =
  ## (sum of a, sum of b, sum of the per-key minimum)
  var (sa, sb, common) = (0, 0, 0)
  for k, n in a:
    sa += n
    common += min(n, b.getOrDefault(k))
  for n in b.values: sb += n
  (sa, sb, common)

func tanimoto*(a, b: CountTable[uint32]): float =
  ## For count fingerprints: Σ min / (Σ a + Σ b − Σ min).
  let (sa, sb, common) = countSums(a, b)
  ratio(common.float, (sa + sb - common).float)

func dice*(a, b: CountTable[uint32]): float =
  ## For count fingerprints: 2 Σ min / (Σ a + Σ b).
  let (sa, sb, common) = countSums(a, b)
  ratio(2.0 * common.float, (sa + sb).float)

func bulkTanimoto*[T: BitVect | CountTable[uint32]](query: T; fps: openArray[T]): seq[float] =
  ## Tanimoto similarity of `query` to each of `fps`.
  result = newSeqOfCap[float](fps.len)
  for fp in fps: result.add tanimoto(query, fp)
