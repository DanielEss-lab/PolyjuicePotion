import openbabel
import pybel
import glob
import math


class Filter:
    def __init__(self, mol):
        self.mol = mol
        for atom in self.mol:
            if atom.OBAtom.IsMetal():
                self.metalIndex = atom.idx - 1
                print(f"Index: {self.metalIndex}; Atom: {self.mol.atoms[self.metalIndex]}")
                print(atom.coords[0])
                break

        # for atom in self.mol:
        #     a = atom.OBAtom
        #     for bond in openbabel.OBAtomBondIter(a):
        #         print(bond.GetId())
        self.calc_distance()

    def calc_distance(self):
        for bond in openbabel.OBAtomBondIter(self.mol.atoms[self.metalIndex].OBAtom):
            print(f"OpenBabel length: {bond.GetLength()}")
