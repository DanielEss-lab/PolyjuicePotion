import openbabel
import pybel
import glob

METAL = 29
BONDED = 7

i = 0
for file in glob.glob("*.xyz"):
    for molecule in pybel.readfile("xyz", file):
        for atom in molecule:
            if atom.OBAtom.IsMetal():
                if atom.atomicnum == METAL:
                    for bond in openbabel.OBAtomBondIter(atom.OBAtom):
                        bonded_to = molecule.atoms[bond.GetNbrAtomIdx(atom.OBAtom) - 1]
                        if bond.GetBondOrder() > 1 and bonded_to.atomicnum == BONDED:
                            print(f"Atom {atom.type}", end=": ")
                            molecule.write("xyz", f"BondLength/{atom.type}{i}.xyz", True)
                            i += 1
                            break
