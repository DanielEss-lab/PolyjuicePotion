# This is a sample Python script.
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import openbabel
import pybel


mol = next(pybel.readfile("xyz", "TestFile.xyz"))
for atom in mol:
    if atom.OBAtom.IsMetal():
        a = atom.OBAtom
        for bond in openbabel.OBAtomBondIter(a):
            print(f"Before bonded to: {mol.atoms[bond.GetNbrAtomIdx(a) - 1]}")
            print(f"Before ID: {bond.GetId()}")
            bond.SetId(1)
            print(f"After bonded to: {mol.atoms[bond.GetNbrAtomIdx(a) - 1]}")
            print(f"After ID: {bond.GetId()}")


mol.write("xyz", "AfterIDChange.xyz", True)
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
