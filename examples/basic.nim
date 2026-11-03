import chem

let aspirin = parseSmiles("O=C(O)c1ccccc1OC(C)=O")
echo "atoms:     ", aspirin.numAtoms, " heavy, ", aspirin.numAtoms(includeH = true), " total"
echo "formula:   ", aspirin.molecularFormula
echo "MW:        ", aspirin.molWeight
echo "exact:     ", aspirin.exactMass
echo "rings:     ", aspirin.numRings
echo "rot bonds: ", aspirin.numRotatableBonds
echo "canonical: ", aspirin.toSmiles
echo "as input:  ", aspirin.toSmiles(canonical = false)

# Build a molecule programmatically: acetonitrile
var m: Mol
let c1 = m.addAtom(initAtom(6))
let c2 = m.addAtom(initAtom(6))
let n = m.addAtom(initAtom(7))
m.addBond(c1, c2)
m.addBond(c2, n, boTriple)
echo m.toSmiles, " ", m.molecularFormula

try:
  discard parseSmiles("C1CC(")
except SmilesError as e:
  echo "error at ", e.pos, ": ", e.msg
