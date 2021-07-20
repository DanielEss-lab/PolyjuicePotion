import time

import openbabel
import pybel
import glob
import re as reeeeeeeeee
import MethylSub


class MonoFinder:
    def __init__(self, mol):
        self.multiBond = False
        self.metalHit = 0
        self.numMetalBonds = int(reeeeeeeeee.findall("MND = (\d+)", str(mol))[0])
        self.mol = mol

        self.find_nearest_atoms()
        self.set_bond_id()

    def nitrogen_multi_bond(self, atom):
        num_bonds = 0
        for _ in openbabel.OBAtomBondIter(atom.OBAtom):
            num_bonds += 1

        if num_bonds < 3:
            return True
        return False

    def find_ligands(self, atom):
        a = atom.OBAtom
        for bond in openbabel.OBAtomBondIter(a):
            if bond.GetId() == 0:
                bond.SetId(1)
                if self.mol.atoms[bond.GetNbrAtomIdx(a) - 1].atomicnum == 1:
                    continue
                elif self.mol.atoms[bond.GetNbrAtomIdx(a) - 1].OBAtom.IsMetal():
                    self.metalHit += 1
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
                    if self.nitrogen_multi_bond(ligand_start):
                        self.metalHit += 1  # If the start of the ligand is a nitrogen with a double or triple bond
                        # the metal, then add one to the metal count to ignore this ligand
                self.find_ligands(ligand_start)
                if self.metalHit == 0:
                    copy_molecule = pybel.Molecule(openbabel.OBMol(self.mol.OBMol))
                    sub = MethylSub.MethylSub(copy_molecule, bond_iter, mol_num)
                    sub.delete_ligand()
                    return  # This causes the first iteration to be the only iteration.  If we want to do all the monodentates, then remove this.
                else:
                    self.metalHit = 0

            bond_iter += 1

    def find_nearest_atoms(self):
        metalBondedList = []
        for atom in self.mol:
            if atom.OBAtom.IsMetal():
                for otherAtom in self.mol:
                    if not otherAtom.OBAtom.IsMetal():
                        metalBondedList.append((otherAtom, atom.OBAtom.GetDistance(otherAtom.OBAtom)))

        newMetalBondedList = sorted(metalBondedList, key=lambda x: x[1])

        for atom in self.mol:
            if atom.OBAtom.IsMetal():
                for j in range(self.numMetalBonds):
                    if self.mol.OBMol.GetBond(atom.OBAtom, newMetalBondedList[j][0].OBAtom) is None:
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
                                self.mol.OBMol.AddBond(bond)
                        else:
                            bond = openbabel.OBBond()
                            bond.SetBegin(atom.OBAtom)
                            bond.SetEnd(newMetalBondedList[j][0].OBAtom)
                            self.mol.OBMol.AddBond(bond)

    def set_bond_id(self):
        for atom in self.mol:
            a = atom.OBAtom
            for bond in openbabel.OBAtomBondIter(a):
                bond.SetId(0)


start_time = time.time()
mol_num = 0
for file in glob.glob("*.xyz"):
    for molecule in pybel.readfile("xyz", file):
        mol_num += 1
        finder = MonoFinder(molecule)
        for atomicBoi in molecule:
            if atomicBoi.OBAtom.IsMetal():
                finder.start(atomicBoi)
end_time = time.time()
elapsed_time = end_time - start_time
minutes = elapsed_time // 60
seconds = elapsed_time % 60
print(f"Ran through {mol_num} molecules in {minutes} minutes and {seconds} seconds.")
