# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 16:38:29 2016

@author: Jonathan Lym
"""
import numpy as np
import warnings
from datetime import datetime, timedelta
from platform import system
from ase.visualize import view
from ase.calculators.vasp import Vasp
from pprint import pprint
import copy

def get_distance(position1, position2, vector = False):
    """
    Provides the distance between two numpy arrays. If vector is True, a numpy array with the
    same dimensions will be returned. Otherwise, a scalar will be returned.
    """
    if vector:
        return position1 - position2
    else:
        return np.sqrt(np.sum([(x1-x2)**2 for x1, x2 in zip(position1, position2)]))

def run_testRun(atoms):
    """
    If running on Linux, prints the ASE Atoms object.
    """
    print('Test Run')
    try:
    	view(atoms)
    except:
    	pass
    pprint(atoms)
        
def print_magmom(atoms):
    """
    Prints the magnetic moment associated with each atom.
    """
    calc = Vasp(istart = 1)
    atoms.set_calculator(calc)
    print('Atom Index\tAtom Type\tMagmom')
    for atom in atoms:
         print('{}}\t{}}\t{}}'.format(atom.index, atom.symbol, atom.magmom))

def print_run_time(out_file):
    """
    For an output file, determines the time taken for each ionic step to complete.
    """
    times = []
    diff_times = []
    i = 0
    with open(out_file, 'r') as out_ptr:
        for line in out_ptr:
            if 'LBFGS' in line:
                times.append(datetime.strptime(line[12:20], '%H:%M:%S'))
                i += 1
    print('Number of steps: {}'.format(i+1))
    diff_times_sec = np.zeros(len(times)-1)
    for i in range(len(times)-1):
        diff_times.append(times[i+1] - times[i])
        if diff_times[i].days < 0:
            diff_times[i] += timedelta(days = 1)
        diff_times_sec[i] = diff_times[i].total_seconds()

    mean_sec = np.mean(diff_times_sec)
    if mean_sec > 3600.:
        print('Mean time per step: {} hours'.format(mean_sec/3600.))
    elif mean_sec > 60.:
        print('Mean time per step: {} minutes'.format(mean_sec/60.))
    else:
        print('Mean time per step: {} seconds'.format(mean_sec))
    tot_sec = np.sum(diff_times_sec)
    if tot_sec > 86400.:
        print('Total time: {} days'.format(tot_sec/86400.))   
    if tot_sec > 3600.:
        print('Total time: {} hours'.format(tot_sec/3600.))
    elif tot_sec > 60.:
        print('Total time: {} minutes'.format(tot_sec/60.))
    else:
        print('Total time: {} seconds'.format(tot_sec))
        
def switch_atoms(atoms, i, j):
    """
    Switches the positions of two atoms. Useful for starting NEB interpolations.
    """
    atoms_copy = atoms.copy()
    atoms[i].position = atoms_copy[j].position
    atoms[j].position = atoms_copy[i].position
    return atoms

def check_bond_lengths(atoms, change_len = False, bond_warning = False, bond_len = 0.9):
    """
    This script goes through an ASE Atoms Object and checks whether any atoms
    are too close. If it finds such a case, it will print a message and change
    the bond length to the critical value
    """
    print("Checking Atoms object: {} for small bond lengths.").format(atoms.get_chemical_formula())
    for atom1 in atoms:
        i1 = atom1.index
        if atom1.symbol != 'H':
            fix = 0
        else:
            fix = 1

        for atom2 in atoms:
            i2 = atom2.index
            if i1 != i2:
                d = atoms.get_distance(i1, i2)
                if d < bond_len:
                    print("Atoms {}({}) and {}({}) have a bond length of {}".format(atom1.symbol, atom1.index, atom2.symbol, atom2.index, d))
                    bond_warning = True
                    if change_len:
                        print("Setting distance to default value of {}".format(bond_len))
                        atoms.set_distance(i1, i2, bond_len, fix = fix)
    if bond_warning == False:
        print("All distances are farther than {}".format(bond_len))

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
