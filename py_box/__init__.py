# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 21:10:42 2016

@author: Jon Lym
"""

def any_alpha(string):
    """Returns True if any alphabetic characters are in the string. False otherwise"""
    for character in string:
        if character.isalpha():
            return True
    else:
        return False

def get_unique_list(data):
    """
    Given a list, returns a unique list.
    """
    keys = {}
    for item in data:
        keys[item] = 1
    return keys.keys()

def switch_atoms(atoms, i, j):
    """
    Switches the positions of two atoms. Useful for starting NEB interpolations.
    """
    atoms_copy = atoms.copy()
    atoms[i].position = atoms_copy[j].position
    atoms[j].position = atoms_copy[i].position
    return atoms

def check_bond_lengths(atoms_obj, change_len = False, bond_warning = False, bond_len = 0.9):
    """
    This script goes through an ASE Atoms Object and checks whether any atoms
    are too close. If it finds such a case, it will print a message and change
    the bond length to the critical value
    """
    print "Checking Atoms object: %s for small bond lengths." % atoms_obj.get_chemical_formula()
    for atom1 in atoms_obj:
        i1 = atom1.index
        if atom1.symbol != 'H':
            fix = 0
        else:
            fix = 1

        for atom2 in atoms_obj:
            i2 = atom2.index
            if i1 != i2:
                d = atoms_obj.get_distance(i1, i2)
                if d < bond_len:
                    print "Atoms %s(%d) and %s(%d) have a bond length of %f" % (atom1.symbol, atom1.index,
                                                                                atom2.symbol, atom2.index,
                                                                                d)
                    bond_warning = True
                    if change_len:
                        print "Setting distance to default value of %f" % bond_len
                        atoms_obj.set_distance(i1, i2, bond_len, fix = fix)
    if bond_warning == False:
        print "All distances are farther than %f" % bond_len
