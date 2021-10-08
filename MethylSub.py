"""
@author Zack Meyer
"""

import os
import pybel
import openbabel
import addH
import LigandChargeFinder
import ChangeCharge


class MethylSub:
    def __init__(self, mol, bond_iter, num_atom):
        self.mol = mol
        self.bond_iter = bond_iter
        self.num_atom = num_atom
        self.metal_ind = 0
        self.new_carbon_ind = 0
        for atom in self.mol:
            atom.OBAtom.SetId(0)
            if atom.OBAtom.IsMetal():
                self.type = atom.type
        self.set_bond_id()
        self.new_charge = self.modify_charge()

    def modify_charge(self):
        for atom in self.mol:
            if atom.OBAtom.IsMetal():
                counter = 0
                for bond in openbabel.OBAtomBondIter(atom.OBAtom):
                    if self.bond_iter == counter:
                        ligand_atom = self.mol.atoms[bond.GetNbrAtomIdx(atom.OBAtom) - 1]
                        charge_finder = LigandChargeFinder.LigandChargeFinder(self.mol)
                        new_charge = charge_finder.change_charge(ligand_atom)

                        # This chunk of code sorts the files into folders based on what their charges are - kept here for future
                        # testing needs
                        # if new_charge == -1:
                        #     if not os.path.exists('Charges'):
                        #         os.makedirs('Charges')
                        #     if not os.path.exists('Charges/MinusOne'):
                        #         os.makedirs('Charges/MinusOne')
                        #     self.mol.write("xyz", f"Charges/MinusOne/{self.type}{self.num_atom}-{self.bond_iter}-Before.xyz", True)
                        # elif new_charge == 1:
                        #     if not os.path.exists('Charges'):
                        #         os.makedirs('Charges')
                        #     if not os.path.exists('Charges/PlusOne'):
                        #         os.makedirs('Charges/PlusOne')
                        #     self.mol.write("xyz", f"Charges/PlusOne/{self.type}{self.num_atom}-{self.bond_iter}-Before.xyz",
                        #                    True)
                        # elif new_charge == 0:
                        #     if not os.path.exists('Charges'):
                        #         os.makedirs('Charges')
                        #     if not os.path.exists('Charges/Zero'):
                        #         os.makedirs('Charges/Zero')
                        #     self.mol.write("xyz", f"Charges/Zero/{self.type}{self.num_atom}-{self.bond_iter}-Before.xyz",
                        #                    True)
                        # else:
                        #     if not os.path.exists('Charges'):
                        #         os.makedirs('Charges')
                        #     if not os.path.exists('Charges/Other'):
                        #         os.makedirs('Charges/Other')
                        #     self.mol.write("xyz", f"Charges/Other/{self.type}{self.num_atom}-{self.bond_iter}-Before.xyz",
                        #                    True)
                        return new_charge

                    counter += 1

    def set_bond_id(self):
        for atom in self.mol:
            a = atom.OBAtom
            for bond in openbabel.OBAtomBondIter(a):
                bond.SetId(0)

    def find_ligands(self, atom):
        a = atom.OBAtom
        a.SetId(2)
        for bond in openbabel.OBAtomBondIter(a):
            if bond.GetId() == 0:
                bond.SetId(1)
                self.find_ligands(self.mol.atoms[bond.GetNbrAtomIdx(a) - 1])

    def delete_ligand(self):
        for atom in self.mol:
            if atom.OBAtom.IsMetal():
                counter = 0
                for bond in openbabel.OBAtomBondIter(atom.OBAtom):
                    bond.SetId(1)
                    if self.bond_iter == counter:
                        self.find_ligands(self.mol.atoms[bond.GetNbrAtomIdx(atom.OBAtom) - 1])
                        # Sets the first atom in the ligand's ID to one so it doesn't get deleted
                        self.mol.atoms[bond.GetNbrAtomIdx(atom.OBAtom) - 1].OBAtom.SetId(1)
                        # Change the atom to a carbon
                        self.mol.atoms[bond.GetNbrAtomIdx(atom.OBAtom) - 1].OBAtom.SetAtomicNum(6)

                        # Set the length to 2.1 so it's close enough to be counted as bonded to the metal
                        bond.SetLength(2.1)

                    counter += 1

        for atom in self.mol:
            if atom.OBAtom.GetId() == 2:
                self.mol.OBMol.DeleteAtom(atom.OBAtom)
                # print("Atom removed")

        # This is after the ligand gets deleted because when they get deleted the indexes of the atoms change
        for atom in self.mol:
            if atom.OBAtom.IsMetal():
                counter = 0
                for bond in openbabel.OBAtomBondIter(atom.OBAtom):
                    if self.bond_iter == counter:
                        # Grab the carbon for the addH function
                        carbon = self.mol.atoms[bond.GetNbrAtomIdx(atom.OBAtom) - 1].OBAtom
                        # Get the index of the carbon needing the hydrogens
                        self.new_carbon_ind = bond.GetNbrAtomIdx(atom.OBAtom)
                        # Get the index of the metal the carbon is bonded to
                        self.metal_ind = bond.GetNbrAtomIdx(carbon)
                    counter += 1

        # This chunk of code sorts the files into folders based on what their charges are - kept here for future
        # testing needs
        # if self.new_charge == -1:
        #     if not os.path.exists('Charges'):
        #         os.makedirs('Charges')
        #     if not os.path.exists('Charges/MinusOne'):
        #         os.makedirs('Charges/MinusOne')
        #     self.mol.write("xyz", f"Charges/MinusOne/{self.type}{self.num_atom}-{self.bond_iter}-After.xyz", True)
        # elif self.new_charge == 1:
        #     if not os.path.exists('Charges'):
        #         os.makedirs('Charges')
        #     if not os.path.exists('Charges/PlusOne'):
        #         os.makedirs('Charges/PlusOne')
        #     self.mol.write("xyz", f"Charges/PlusOne/{self.type}{self.num_atom}-{self.bond_iter}-After.xyz",
        #                    True)
        # elif self.new_charge == 0:
        #     if not os.path.exists('Charges'):
        #         os.makedirs('Charges')
        #     if not os.path.exists('Charges/Zero'):
        #         os.makedirs('Charges/Zero')
        #     self.mol.write("xyz", f"Charges/Zero/{self.type}{self.num_atom}-{self.bond_iter}-After.xyz",
        #                    True)
        # else:
        #     if not os.path.exists('Charges'):
        #         os.makedirs('Charges')
        #     if not os.path.exists('Charges/Other'):
        #         os.makedirs('Charges/Other')
        #     self.mol.write("xyz", f"Charges/Other/{self.type}{self.num_atom}-{self.bond_iter}-After.xyz",
        #                    True)

        if not os.path.exists('DeletedMono'):
            os.makedirs('DeletedMono')
        self.mol.write("xyz", f"DeletedMono/{self.type}{self.num_atom}-{self.bond_iter}.xyz", True)
        charge_changer = ChangeCharge.ChargeChanger(f"DeletedMono/{self.type}{self.num_atom}-{self.bond_iter}.xyz")
        charge_changer.change(self.new_charge)

        hydrogen_adder = addH.AddH(f"DeletedMono/{self.type}{self.num_atom}-{self.bond_iter}.xyz", self.metal_ind, self.new_carbon_ind)
        hydrogen_adder.start()
