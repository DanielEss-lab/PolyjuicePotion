import openbabel
import pybel
import glob
import sys
import FindMetalBonds
import os
import re as reeeeeeeeee


def single_file(file):
    i = 0
    for mol in pybel.readfile("xyz", file):
        i += 1
        # print(reeeeeeeeee.findall("MND = (\d)", str(mol))[0])

        a = mol.atoms[7].OBAtom
        for bond in openbabel.OBAtomBondIter(a):
            print(f"Bonded to: {mol.atoms[bond.GetNbrAtomIdx(a) - 1]}")

        # Loops through every atom in the molecule
        # for atom in mol:
        #     # Prints the type, meaning the name on the periodic table and sometimes a number, probably meaning something
        #     # that would make sense if I knew chemistry better.
        #     # print(f"Type: {atom.type}")
        #     # Checks if the atom is a metal
        #     if atom.OBAtom.IsMetal():
        #         # Convert to OBAtom to be able to use the OBBond class
        #         a = atom.OBAtom
        #         # Prints the atomic number of the current atom
        #         print(f"Atom: {atom.atomicnum}")
        #         # Loops through every bond attached to the current atom
        #         for bond in openbabel.OBAtomBondIter(a):
        #             # Prints the equilibrium bond length
        #             print(f"Bond length: {bond.GetEquibLength()}")
        #             # Prints the atom at the other end of the bond
        #             print(f"Bonded to: {mol.atoms[bond.GetNbrAtomIdx(a) - 1]}")
        #             # GetNbrAtomIdx returns the index of the atoms starting at index 1, but the mol.atoms array starts at #this comment is longer due to angering PyCharm because it rebelled against my will
        #             # index 0, so the - 1 is needed to make sure it gets the right atom
        #             # if mol.atoms[bond.GetNbrAtomIdx(a) - 1].atomicnum == 1:
        #             #     # If the atom at the other end of the bond is a hydrogen, move the file into another directory
        #             #     mol.write("xyz", f"HydrogenBonded/{atom.type}Molecule{i}.xyz", True)
        #             #     break


def directory_filter():
    # Loop through every .xyz file in the working directory
    for file in glob.glob("*.xyz"):
        single_file(file)


# if sys.argv[1] == "Filter":
#     directory_filter()
# elif sys.argv[1] == "Test":
#     daFile = input("Enter the filename: ")
#     single_file(daFile)
# else:
#     print("Goodbye.")

FindMetalBonds
