import openbabel
import pybel
import glob
import TestSwitch
import re as reeeeeeeeee


METAL = 74
BONDED = 8
BOND_LENGTH = 1.8


def find_nearest_atoms(mol):
    numMetalBonds = int(reeeeeeeeee.findall("MND = (\d+)", str(mol))[0])
    metalBondedList = []
    for atom in mol:
        if atom.OBAtom.IsMetal():
            for otherAtom in mol:
                if not otherAtom.OBAtom.IsMetal():
                    metalBondedList.append((otherAtom, atom.OBAtom.GetDistance(otherAtom.OBAtom)))

    newMetalBondedList = sorted(metalBondedList, key=lambda x: x[1])

    for atom in mol:
        if atom.OBAtom.IsMetal():
            for j in range(numMetalBonds):
                if mol.OBMol.GetBond(atom.OBAtom, newMetalBondedList[j][0].OBAtom) is None:
                    if newMetalBondedList[j][0].atomicnum == 1:
                        i = 0
                        for _ in openbabel.OBAtomBondIter(newMetalBondedList[j][0].OBAtom):
                            i += 1
                        if i > 0:
                            break
                        else:
                            bond = openbabel.OBBond()
                            bond.SetBegin(atom.OBAtom)
                            bond.SetEnd(newMetalBondedList[j][0].OBAtom)
                            mol.OBMol.AddBond(bond)
                    else:
                        bond = openbabel.OBBond()
                        bond.SetBegin(atom.OBAtom)
                        bond.SetEnd(newMetalBondedList[j][0].OBAtom)
                        mol.OBMol.AddBond(bond)
                        # print(pybel.Atom(bond.GetBeginAtom()))
                        # print(pybel.Atom(bond.GetEndAtom()))


i = 0
for file in glob.glob("*.xyz"):
    for molecule in pybel.readfile("xyz", file):
        find_nearest_atoms(molecule)
        for atom in molecule:
            if atom.OBAtom.IsMetal():
                # if atom.atomicnum == METAL:
                for bond_atom in openbabel.OBAtomBondIter(atom.OBAtom):
                    bonded_to = molecule.atoms[bond_atom.GetNbrAtomIdx(atom.OBAtom) - 1]
                    # print(bonded_to)
                    if bonded_to.atomicnum == BONDED \
                            and bond_atom.GetLength() < BOND_LENGTH:
                        # print(f"Atom {atom.type}", end=": ")
                        TestSwitch.switch(atom.atomicnum, bonded_to.atomicnum)
                        molecule.write("xyz", f"BondLength/{atom.type}{i}.xyz", True)
                        i += 1
                        break
