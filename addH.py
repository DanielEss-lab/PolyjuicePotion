# Created by Nathan Morgan June 2021
# Usage: python addH.py coordinates.xyz
# Adds hydrogen to carbon to make a methyl group
# Make sure you use python 3 (module load python/3.8) 
import csv
import numpy as np
import os

# Filename = sys.argv[1]  # Get the name of the structure to edit from the commandline


class AddH:
    def __init__(self, filename, metal_ind, carbon_ind):
        self.filename = filename
        self.metal_ind = metal_ind
        self.carbon_ind = carbon_ind

    def readxyz(self):
        atoms = np.empty(0)  # Initialize array for atom elements
        coord = np.zeros((0, 3))  # and cartesian coordinates

        with open(self.filename, 'r', newline='') as f:
            reader = csv.reader(f, delimiter=' ')
            for row in reader:
                row = [i for i in row if i]  # Remove whitespace and put elements in a list
                if len(row) == 4:  # Lines like "C 1.23 4.56 7.89" have 4 elements
                    # print(f"Row: {row}")
                    atoms = np.append(atoms, [row[0]])
                    coord = np.append(coord, [row[1:4]], axis=0)
        return atoms, coord.astype(np.float)  # Return the coordinates as numbers not strings

    def writexyz(self, coord1, coord2, coord3):
        orig = open(self.filename, 'r', newline='')
        if not os.path.exists('HydrogensAdded'):
            os.makedirs('HydrogensAdded')
        new = open(f'HydrogensAdded/H-{self.filename[12:]}', 'w')
        reader = csv.reader(orig, delimiter=' ')
        for row in reader:
            row = [i for i in row if i]
            if len(row) == 1:  # The first line contains only the number of atoms
                new.write(str(int(row[0]) + 3))  # Add three to the number of atoms because we add three hydrogens
                new.write('\n')
            else:
                string = ' '.join([str(elem) for elem in row])  # I'm not sure why I needed to make this a string
                new.write(string)  # but I guess I did
                if len(row) > 4:    # If the row has more than just the atom and its coords, add the carbon index
                    new.write(f" | CI: {self.carbon_ind}")
                new.write('\n')
        string1 = ' '.join([str(elem) for elem in coord1])  # These coordinates need to be a string so I can
        string2 = ' '.join([str(elem) for elem in coord2])  # concatenate them
        string3 = ' '.join([str(elem) for elem in coord3])
        new.write('H ' + string1)
        new.write('\n')
        new.write('H ' + string2)
        new.write('\n')
        new.write('H ' + string3)
        orig.close()
        new.close()
        return

    def getvec(self, coord, atom1, atom2):  # Returns a vector between two atoms in a list of coordinates
        # print(f"Atom Index: {atom1}\nCoords: {coord[atom1]}")
        # print(f"Atom Index: {atom2}\nCoords: {coord[atom2]}")
        [x1, y1, z1] = coord[atom1 - 1]  # I let the atoms be 1-indexed so minus 1 to convert to
        [x2, y2, z2] = coord[atom2 - 1]  # zero indexed python
        return np.array([x2 - x1, y2 - y1, z2 - z1])

    def getangle(self, vec1, vec2):  # This function returns the angle between two vectors. It isn't actually used.
        c = np.dot(vec1, vec2) / np.linalg.norm(vec1) / np.linalg.norm(vec2)
        return np.arccos(np.clip(c, -1, 1)) * 180 / np.pi

                                                # H
                                            #	 /
    def bit_to_add(self, vec):  # In methane H--C-H, the three hydrogens stick out parallel to the first hydrogen by 0.357 angstroms
        length = np.linalg.norm(vec)           # \
        return 0.357 * vec / length             # H

    def getperp(self, vec):
        try1 = np.array(
            [vec[1], -vec[0], 0])  # These three vectors are all perpendicular to the given vector by construction
        try2 = np.array([0, vec[2], -vec[1]])  # (forces the dot product to be 0)
        try3 = np.array(
            [vec[2], 0,
             -vec[0]])  # Give it three trys because a vector like <0,0,1> would give a perpendicular vector of
        if np.sum(try1) != 0.:  # <0,0,0> using try1's method
            return try1  # return the first perpedicular vector that isn't the zero vector
        if np.sum(try2) != 0.:
            return try2
        if np.sum(try3) != 0.:
            return try3

    def perpperp(self, vec1, vec2):  # Returns a vector perpendicular to two given vector using the cross product
        [a1, a2, a3] = vec1
        [b1, b2, b3] = vec2
        return np.array([(a2 * b3 - a3 * b2), (a3 * b1 - a1 * b3), (a1 * b2 - a2 * b1)])

    def go_up(self, vec):
        length = np.linalg.norm(vec)  # H			#	H
        return 1.01 * vec / length  # /|			#       |

    def go_down(self, vec):  # / |1.01 Angstrom	#       |
        length = np.linalg.norm(vec)  # H-----C  |			#     H>C
        return 0.504 * vec / length  # \\   |			#      / \

    def go_over(self, vec):  # HH  |0.504 Angstrom	#     /   \
        length = np.linalg.norm(vec)  # H     H
        return 0.878 * vec / length  # --- 0.878 Angstrom

    def start(self):
        atoms, coord = self.readxyz()  # Read in the coordinates
        # print(atoms)				# Uncomment to help debug
        # print(coord)
        vec = self.getvec(coord, self.metal_ind, self.carbon_ind)
        # This assumes the metal is the first atom in the list and the carbon to add hydrogen to is the second
        # Could be adjusted by changing the atom numbers

        coord1 = np.round(coord[self.carbon_ind - 1] + self.bit_to_add(vec) + self.go_up(self.getperp(vec)), 5)
        # coord[self.carbon_ind - 1] is the position of the carbon and the additions shift the hydrogen into place
        coord2 = coord[self.carbon_ind - 1] + self.bit_to_add(vec) - self.go_down(self.getperp(vec)) + self.go_over(
            self.perpperp(vec, self.getperp(vec)))  # Could also round these ones if desired
        coord3 = coord[self.carbon_ind - 1] + self.bit_to_add(vec) - self.go_down(self.getperp(vec)) - self.go_over(self.perpperp(vec, self.getperp(vec)))

        self.writexyz(coord1, coord2, coord3)  # Write the output file with the three hydrogens appeneded
