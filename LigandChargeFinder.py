import openbabel
import pybel


class LigandChargeFinder:
    def __init__(self, molecule):
        self.molecule = molecule

    def bondcounter(self, a):
        bond_count = 0
        for _ in openbabel.OBAtomBondIter(a):
            bond_count += 1
        return bond_count

    def cofinder(self, carbon):
        c = carbon.OBAtom
        carbon_bond_count = self.bondcounter(c)
        for bond in openbabel.OBAtomBondIter(c):
            atomic_num = self.molecule.atoms[bond.GetNbrAtomIdx(c) - 1].atomicnum
            if atomic_num != 8 and atomic_num != 16:
                continue
            if carbon_bond_count == 2:
                oxygen = self.molecule.atoms[bond.GetNbrAtomIdx(c) - 1].OBAtom
                oxygen_bond_count = self.bondcounter(oxygen)
                if oxygen_bond_count == 1:
                    return True
                return False
            return False
        return False

    def cnfinder(self, carbon):
        c = carbon.OBAtom
        carbon_bond_count = self.bondcounter(c)
        for bond in openbabel.OBAtomBondIter(c):
            atom_x = self.molecule.atoms[bond.GetNbrAtomIdx(c) - 1].atomicnum
            if atom_x != 7:
                continue
            if carbon_bond_count == 2:
                nitrogen = self.molecule.atoms[bond.GetNbrAtomIdx(c) - 1].OBAtom
                nitrogen_bond_count = self.bondcounter(nitrogen)
                if nitrogen_bond_count > 1:
                    return True
                return False
            return False
        return False

    def wateresque_finder(self, oxygen):
        o = oxygen.OBAtom
        oxygen_bond_count = self.bondcounter(o)
        if oxygen_bond_count > 2 and (oxygen.atomicnum == 8 or oxygen.atomicnum == 16):
            return True
        return False

    def ammonia_finder(self, nitrogen):
        n = nitrogen.OBAtom
        nitrogen_bond_count = self.bondcounter(n)
        if nitrogen_bond_count > 3 and (nitrogen.atomicnum == 7 or nitrogen.atomicnum == 15):
            return True
        return False

    def oxygen_double_bond(self, oxygen):
        o = oxygen.OBAtom
        bonds = self.bondcounter(o)
        if bonds == 1 and oxygen.atomicnum == 8:
            return True
        return False

    def carbon_double_bond(self, carbon):
        c = carbon.OBAtom
        bonds = self.bondcounter(c)
        if bonds == 3 and carbon.atomicnum == 6:
            return True
        if bonds == 2 and carbon.atomicnum == 6:
            if self.cofinder(carbon) or self.cnfinder(carbon):
                return False
            return True
        return False

    def help_aromatic_ring(self, atom, iteration):
        num_bonds = self.bondcounter(atom.OBAtom)
        if (atom.atomicnum == 6 and num_bonds != 3) or (atom.atomicnum == 7 and num_bonds > 3):
            return False
        if iteration > 6:
            return False
        a = atom.OBAtom
        for bond in openbabel.OBAtomBondIter(a):
            other_end = self.molecule.atoms[bond.GetNbrAtomIdx(a) - 1]
            if other_end.OBAtom.GetId() == 3 and bond.GetId() == 1:
                return True
            if bond.GetId() == 1:
                bond.SetId(2)
                return self.help_aromatic_ring(other_end, iteration + 1)

        return False

    def find_aromatic_ring(self, atom):
        if atom.atomicnum != 6 and atom.atomicnum != 7:     # The aromatic ring should start with either a nitrogen or a carbon
            return False
        num_bonds = self.bondcounter(atom.OBAtom)
        if (atom.atomicnum == 6 and num_bonds != 3) or (atom.atomicnum == 7 and num_bonds > 3):
            return False
        for bond in openbabel.OBMolBondIter(self.molecule.OBMol):
            bond.SetId(1)
        for element in self.molecule:
            element.OBAtom.SetId(1)
        atom.OBAtom.SetId(3)
        a = atom.OBAtom
        for bond in openbabel.OBAtomBondIter(a):
            if self.molecule.atoms[bond.GetNbrAtomIdx(a) - 1].OBAtom.IsMetal():
                bond.SetId(2)
            if bond.GetId() == 1:
                bond.SetId(2)
                is_ring = self.help_aromatic_ring(self.molecule.atoms[bond.GetNbrAtomIdx(a) - 1], 1)

        return is_ring

    def charge_change(self, atom):
        neutral_ligand = self.cnfinder(atom) or self.cofinder(atom) or self.wateresque_finder(atom) or \
                         self.ammonia_finder(atom)
        double_bond = self.oxygen_double_bond(atom) or self.carbon_double_bond(atom)
        triple_bond = False
        if neutral_ligand:
            return -1
        elif double_bond:
            return 1
        elif triple_bond:
            return 2
        else:                   # Negative ligands
            return 0
