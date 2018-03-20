# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 16:38:29 2016

@author: Jonathan Lym
"""

from py_box.ase.atom import Atom
from py_box.ase.atoms import Atoms
import numpy as np
import warnings
from datetime import datetime, timedelta
from platform import system
from pprint import pprint
from ase.visualize import view
from ase.calculators.vasp import Vasp
import copy

__all__ = ['Atoms', 'Atom']

def get_distance(position1, position2, vector = False):
    if vector:
        return position1 - position2
    else:
        return np.sqrt(np.sum([(x1-x2)**2 for x1, x2 in zip(position1, position2)]))

def run_testRun(atoms_obj):
    print 'Test Run'
    os_name = system()
#    if type(atoms_obj) is list:        
#        for i in range(len(atoms_obj)):    
#            pprint(vars(atoms_obj[i]))
#    else:
#        pprint(vars(atoms_obj))
    try:
        view(atoms_obj)
    except:
        pass
        
def print_magmom(atoms_obj):
    """
    Prints the magnetic moment associated with each atom.
    """
    calc = Vasp(istart = 1)
    atoms_obj.set_calculator(calc)
    print 'Atom Index\tAtom Type\tMagmom'
    for atom in atoms_obj:
         print '%d\t%s\t%f' % (atom.index, atom.symbol, atom.magmom)

def print_run_time(out_file):
    times = []
    diff_times = []
    i = 0
    with open(out_file, 'r') as out_ptr:
        for line in out_ptr:
            if 'LBFGS' in line:
                times.append(datetime.strptime(line[12:20], '%H:%M:%S'))
                i += 1
    print 'Number of steps: %d' % (i+1)
    diff_times_sec = np.zeros(len(times)-1)
    for i in range(len(times)-1):
        diff_times.append(times[i+1] - times[i])
        if diff_times[i].days < 0:
            diff_times[i] += timedelta(days = 1)
        diff_times_sec[i] = diff_times[i].total_seconds()

    mean_sec = np.mean(diff_times_sec)
    if mean_sec > 3600.:
        print 'Mean time per step: %f hours' % (mean_sec/3600.)
    elif mean_sec > 60.:
        print 'Mean time per step: %f minutes' % (mean_sec/60.)
    else:
        print 'Mean time per step: %f seconds' % mean_sec
        
    tot_sec = np.sum(diff_times_sec)
    if tot_sec > 86400.:
        print 'Total time: %f days' % (tot_sec/86400.)        
    if tot_sec > 3600.:
        print 'Total time: %f hours' % (tot_sec/3600.)
    elif tot_sec > 60.:
        print 'Total time: %f minutes' % (tot_sec/60.)
    else:
        print 'Total time: %f seconds' % tot_sec
        
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

DFT_E_gas = {
    'H2': -6.759196,
    'H2O': -14.219543,
    'FAL': -79.56456956,
    'MF': -73.257932

}

# def insert_atom(atoms, atom, index):
#     """
#     Inserts an atom into the atoms object and shifts all the indices.
#     """
#     offset = 0
#     out_atoms = atoms.copy().append(atom)
#     for i in range(len(atoms)+1):
#         if i == index:
#             out_atoms[i] = copy.copy(atom)
#             offset += 1
#         else:
#             out_atoms[i+offset] = copy.copy(atoms[i])
#     return out_atoms
