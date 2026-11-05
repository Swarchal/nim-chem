## Periodic table data: symbols, masses and default valences.
##
## Average masses are IUPAC standard atomic weights (conventional values);
## monoisotopic masses are those of the most abundant isotope.

import std/tables
import ./isotopes

type ElementInfo* = object
  symbol*: string
  avgMass*: float ## standard atomic weight
  monoMass*: float ## mass of the most abundant isotope
  valences*: seq[int] ## allowed valences (empty = no default valence model)

func e(sym: string, avg, mono: float, vals: seq[int] = @[]): ElementInfo =
  ElementInfo(symbol: sym, avgMass: avg, monoMass: mono, valences: vals)

const elementTable*: array[0 .. 118, ElementInfo] = [
  e("*", 0.0, 0.0),
  e("H", 1.008, 1.00782503207, @[1]),
  e("He", 4.002602, 4.00260325415, @[0]),
  e("Li", 6.94, 7.01600455, @[1]),
  e("Be", 9.0121831, 9.0121822, @[2]),
  e("B", 10.81, 11.0093054, @[3]),
  e("C", 12.011, 12.0, @[4]),
  e("N", 14.007, 14.0030740048, @[3, 5]),
  e("O", 15.999, 15.99491461956, @[2]),
  e("F", 18.998403163, 18.99840322, @[1]),
  e("Ne", 20.1797, 19.9924401754, @[0]),
  e("Na", 22.98976928, 22.9897692809, @[1]),
  e("Mg", 24.305, 23.985041700, @[2]),
  e("Al", 26.9815385, 26.98153863, @[3]),
  e("Si", 28.085, 27.9769265325, @[4]),
  e("P", 30.973761998, 30.97376163, @[3, 5]),
  e("S", 32.06, 31.97207100, @[2, 4, 6]),
  e("Cl", 35.45, 34.96885268, @[1, 3, 5, 7]),
  e("Ar", 39.948, 39.9623831225, @[0]),
  e("K", 39.0983, 38.96370668, @[1]),
  e("Ca", 40.078, 39.96259098, @[2]),
  e("Sc", 44.955908, 44.9559119),
  e("Ti", 47.867, 47.9479463),
  e("V", 50.9415, 50.9439595),
  e("Cr", 51.9961, 51.9405075),
  e("Mn", 54.938044, 54.9380451),
  e("Fe", 55.845, 55.9349375),
  e("Co", 58.933194, 58.9331950),
  e("Ni", 58.6934, 57.9353429),
  e("Cu", 63.546, 62.9295975),
  e("Zn", 65.38, 63.9291422),
  e("Ga", 69.723, 68.9255736, @[3]),
  e("Ge", 72.630, 73.9211778, @[4]),
  e("As", 74.921595, 74.9215965, @[3, 5]),
  e("Se", 78.971, 79.9165213, @[2, 4, 6]),
  e("Br", 79.904, 78.9183371, @[1, 3, 5, 7]),
  e("Kr", 83.798, 83.911507, @[0]),
  e("Rb", 85.4678, 84.911789738, @[1]),
  e("Sr", 87.62, 87.9056121, @[2]),
  e("Y", 88.90584, 88.9058483),
  e("Zr", 91.224, 89.9047044),
  e("Nb", 92.90637, 92.9063781),
  e("Mo", 95.95, 97.9054082),
  e("Tc", 98.0, 97.9072160),
  e("Ru", 101.07, 101.9043493),
  e("Rh", 102.90550, 102.905504),
  e("Pd", 106.42, 105.903486),
  e("Ag", 107.8682, 106.905097),
  e("Cd", 112.414, 113.9033585),
  e("In", 114.818, 114.903878, @[3]),
  e("Sn", 118.710, 119.9021947, @[2, 4]),
  e("Sb", 121.760, 120.9038157, @[3, 5]),
  e("Te", 127.60, 129.9062244, @[2, 4, 6]),
  e("I", 126.90447, 126.904473, @[1, 3, 5, 7]),
  e("Xe", 131.293, 131.9041535),
  e("Cs", 132.90545196, 132.905451933, @[1]),
  e("Ba", 137.327, 137.9052472, @[2]),
  e("La", 138.90547, 138.9063533),
  e("Ce", 140.116, 139.9054387),
  e("Pr", 140.90766, 140.9076528),
  e("Nd", 144.242, 141.9077233),
  e("Pm", 145.0, 144.912749),
  e("Sm", 150.36, 151.9197324),
  e("Eu", 151.964, 152.9212303),
  e("Gd", 157.25, 157.9241039),
  e("Tb", 158.92535, 158.9253468),
  e("Dy", 162.500, 163.9291748),
  e("Ho", 164.93033, 164.9303221),
  e("Er", 167.259, 165.9302931),
  e("Tm", 168.93422, 168.9342133),
  e("Yb", 173.045, 173.9388621),
  e("Lu", 174.9668, 174.9407718),
  e("Hf", 178.49, 179.9465500),
  e("Ta", 180.94788, 180.9479958),
  e("W", 183.84, 183.9509312),
  e("Re", 186.207, 186.9557531),
  e("Os", 190.23, 191.9614807),
  e("Ir", 192.217, 192.9629264),
  e("Pt", 195.084, 194.9647911),
  e("Au", 196.966569, 196.9665687),
  e("Hg", 200.592, 201.9706430),
  e("Tl", 204.38, 204.9744275),
  e("Pb", 207.2, 207.9766521),
  e("Bi", 208.98040, 208.9803987),
  e("Po", 209.0, 208.9824304),
  e("At", 210.0, 209.987148),
  e("Rn", 222.0, 222.0175777),
  e("Fr", 223.0, 223.0197359),
  e("Ra", 226.0, 226.0254098),
  e("Ac", 227.0, 227.0277521),
  e("Th", 232.0377, 232.0380553),
  e("Pa", 231.03588, 231.0358840),
  e("U", 238.02891, 238.0507882),
  e("Np", 237.0, 237.0481734),
  e("Pu", 244.0, 244.064204),
  e("Am", 243.0, 243.0613811),
  e("Cm", 247.0, 247.070354),
  e("Bk", 247.0, 247.070307),
  e("Cf", 251.0, 251.079587),
  e("Es", 252.0, 252.082980),
  e("Fm", 257.0, 257.095105),
  e("Md", 258.0, 258.098431),
  e("No", 259.0, 259.10103),
  e("Lr", 266.0, 262.10963),
  e("Rf", 267.0, 267.0),
  e("Db", 268.0, 268.0),
  e("Sg", 269.0, 269.0),
  e("Bh", 270.0, 270.0),
  e("Hs", 269.0, 269.0),
  e("Mt", 278.0, 278.0),
  e("Ds", 281.0, 281.0),
  e("Rg", 282.0, 282.0),
  e("Cn", 285.0, 285.0),
  e("Nh", 286.0, 286.0),
  e("Fl", 289.0, 289.0),
  e("Mc", 290.0, 290.0),
  e("Lv", 293.0, 293.0),
  e("Ts", 294.0, 294.0),
  e("Og", 294.0, 294.0),
]

const maxAtomicNumber* = elementTable.high

func isValidElement*(z: int): bool {.inline.} =
  z in 0 .. maxAtomicNumber

func symbol*(z: int): string =
  elementTable[z].symbol
func avgMass*(z: int): float =
  elementTable[z].avgMass
func monoMass*(z: int): float =
  elementTable[z].monoMass
func defaultValences*(z: int): seq[int] =
  elementTable[z].valences

func buildSymbolTable(): Table[string, int] =
  for z in 0 .. maxAtomicNumber:
    result[elementTable[z].symbol] = z

const symbolTable = buildSymbolTable()

func atomicNumber*(sym: string): int =
  ## Atomic number for an element symbol (case-sensitive, e.g. "Cl"), or -1.
  symbolTable.getOrDefault(sym, -1)

func isotopeMass*(z, massNumber: int): float =
  ## Exact mass of a specific isotope (3,111 known, see `isotopes`). Falls
  ## back to the mass number for isotopes not in the table.
  if z notin isotopeMasses.low .. isotopeMasses.high:
    return float(massNumber)
  let (first, masses) = isotopeMasses[z]
  let k = massNumber - first
  if k in 0 .. masses.high:
    masses[k]
  else:
    float(massNumber)

const
  organicSubset* = ["B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"]
  aromaticSymbols* = ["b", "c", "n", "o", "p", "s", "se", "as", "te"]
    ## element symbols that may be written in lowercase (aromatic) form
  organicAromatic* = ["b", "c", "n", "o", "p", "s"]
    ## lowercase symbols allowed outside of brackets

  organicSubsetZ = {5'u8, 6, 7, 8, 9, 15, 16, 17, 35, 53}
  organicAromaticZ = {5'u8, 6, 7, 8, 15, 16}
  aromaticZ = organicAromaticZ + {33'u8, 34, 52}

func inOrganicSubset*(z: int): bool =
  uint8(z) in organicSubsetZ
func canBeAromatic*(z: int): bool =
  uint8(z) in aromaticZ
func inOrganicAromatic*(z: int): bool =
  uint8(z) in organicAromaticZ
