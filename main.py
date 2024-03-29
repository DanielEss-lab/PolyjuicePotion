"""
@author Zack Meyer
"""

import time

from openbabel import openbabel
from openbabel import pybel
import glob
import re
import MethylSub


def find_nitrogen_multi_bond(atom):
    num_bonds = 0
    for _ in openbabel.OBAtomBondIter(atom.OBAtom):
        num_bonds += 1

    if num_bonds < 3:
        return True
    return False


def find_num_bonds_to_metal(mol):
    mnd_list = re.findall("MND = (\d+)", str(mol));
    return mnd_list[0] if mnd_list != [] else "0"


class MonoFinder:
    def __init__(self, mol):
        self.multi_bond = False
        self.metal_hit = 0
        self.num_metal_bonds = int(find_num_bonds_to_metal(mol))
        self.mol = mol

        self.find_nearest_atoms()
        self.set_bond_id()

    def find_ligands(self, atom):
        a = atom.OBAtom
        for bond in openbabel.OBAtomBondIter(a):
            if bond.GetId() == 0:
                bond.SetId(1)
                if self.mol.atoms[bond.GetNbrAtomIdx(a) - 1].atomicnum == 1:
                    continue
                elif self.mol.atoms[bond.GetNbrAtomIdx(a) - 1].OBAtom.IsMetal():
                    self.metal_hit += 1
                    continue
                else:
                    self.find_ligands(self.mol.atoms[bond.GetNbrAtomIdx(a) - 1])

    def start(self, metal):
        bond_iter = 0
        m = metal.OBAtom
        for bond in openbabel.OBAtomBondIter(m):
            if bond.GetId() == 0:
                bond.SetId(1)
                ligand_start = self.mol.atoms[bond.GetNbrAtomIdx(m) - 1]
                if ligand_start.atomicnum == 7:
                    if find_nitrogen_multi_bond(ligand_start):
                        self.metal_hit += 1  # If the start of the ligand is a nitrogen with a double or triple bond
                        # the metal, then add one to the metal count to ignore this ligand

                if not do_ligand_filter:
                    self.find_ligands(ligand_start)
                elif do_ligand_filter and ligand_start.atomicnum == ligand_filter:
                    self.find_ligands(ligand_start)
                else:
                    bond_iter += 1
                    continue
                if self.metal_hit == 0:
                    copy_molecule = pybel.Molecule(openbabel.OBMol(self.mol.OBMol))
                    sub = MethylSub.MethylSub(copy_molecule, bond_iter, mol_num)
                    sub.delete_ligand()
                    return  # This causes the first iteration to be the only iteration.  If we want to do all the monodentates, then remove this.
                else:
                    self.metal_hit = 0

            bond_iter += 1

    def find_nearest_atoms(self):
        metal_bonded_list = []
        for atom in self.mol:
            if atom.OBAtom.IsMetal():
                for otherAtom in self.mol:
                    if not otherAtom.OBAtom.IsMetal():
                        metal_bonded_list.append((otherAtom, atom.OBAtom.GetDistance(otherAtom.OBAtom)))

        new_metal_bonded_list = sorted(metal_bonded_list, key=lambda x: x[1])

        for atom in self.mol:
            if atom.OBAtom.IsMetal():
                for j in range(self.num_metal_bonds):
                    if self.mol.OBMol.GetBond(atom.OBAtom, new_metal_bonded_list[j][0].OBAtom) is None:
                        if new_metal_bonded_list[j][0].atomicnum == 1:
                            i = 0
                            for _ in openbabel.OBAtomBondIter(new_metal_bonded_list[j][0].OBAtom):
                                i += 1
                            if i > 0:
                                break
                            else:
                                bond = openbabel.OBBond()
                                bond.SetBegin(atom.OBAtom)
                                bond.SetEnd(new_metal_bonded_list[j][0].OBAtom)
                                self.mol.OBMol.AddBond(bond)
                        else:
                            bond = openbabel.OBBond()
                            bond.SetBegin(atom.OBAtom)
                            bond.SetEnd(new_metal_bonded_list[j][0].OBAtom)
                            self.mol.OBMol.AddBond(bond)

    def set_bond_id(self):
        for atom in self.mol:
            a = atom.OBAtom
            for bond in openbabel.OBAtomBondIter(a):
                bond.SetId(0)


filter_input = input("Do you want to only replace ligands for a certain metal? ").lower()
do_metal_filter = filter_input == 'y' or filter_input == 'yes'
metal_filter = 0
if do_metal_filter:
    metal_filter = int(input("Please enter the atomic number of the metal you want: "))

filter_input = input("Do you want to only replace ligands with a certain starting atom? ").lower()
do_ligand_filter = filter_input == 'y' or filter_input == 'yes'
ligand_filter = 0
if do_ligand_filter:
    ligand_filter = int(input("Please enter the atomic number of the ligand's starting atom: "))

start_time = time.time()
mol_num = 0
for file in glob.glob("*.xyz"):
    for molecule in pybel.readfile("xyz", file):
        mol_num += 1
        finder = MonoFinder(molecule)
        for metal in molecule:
            if metal.OBAtom.IsMetal() and not do_metal_filter:
                finder.start(metal)
            elif metal.OBAtom.IsMetal() and metal.atomicnum == metal_filter:
                finder.start(metal)
            # else:
            #     break

end_time = time.time()
elapsed_time = end_time - start_time
minutes = elapsed_time // 60
seconds = elapsed_time % 60
print(f"Ran through {mol_num} molecules in {minutes:.0f} minutes and {seconds:.2f} seconds.")
