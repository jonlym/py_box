# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 16:38:29 2016

@author: Jonathan Lym
"""

import numpy as np
from platform import system
from pprint import pprint
from ase.visualize import view
from ase.calculators.vasp import Vasp

def run_testRun(atoms_obj):
    print 'Test Run'
    os_name = system()
#    if type(atoms_obj) is list:        
#        for i in range(len(atoms_obj)):    
#            pprint(vars(atoms_obj[i]))
#    else:
#        pprint(vars(atoms_obj))
    if os_name.lower() == 'linux':
        view(atoms_obj)
        
def find_coordination_number(atoms_obj, max_len = 2.5):
    """
    Determines the coordination number based on the max_len constraint. Note 
    that ASE has a neighbors class (ase.neighborlist), which may be more useful
    than this script
    """
    n_atoms = len(atoms_obj)
    #Duplicate the atoms object so that it is surrounded to get an accurate CN
    atoms_ext = atoms_obj.copy()
    atoms_ext *= (3, 3, 1)

    offset = 4 * n_atoms
    for i in range(n_atoms):
        print '-'*20
        print 'Finding nearest neighbors to atom %s [%d]' % (atoms_obj[i].symbol, i)
        CN = 0
        for j in range(len(atoms_ext)):    
            if (i + offset) != j:
                bond_len = atoms_ext.get_distance(i + offset, j)
                if bond_len < max_len:
                    k = int(j - n_atoms*np.floor(j / n_atoms))
                    print '\t%s [%d]\t %f' % (atoms_obj[k].symbol, k, bond_len)
                    CN += 1
        print '\t Coordination Number: %d' % CN

def print_magmom(atoms_obj):
    """
    Prints the magnetic moment associated with each atom.
    """
    calc = Vasp(istart = 1)
    atoms_obj.set_calculator(calc)
    print 'Atom Index\tAtom Type\tMagmom'
    for atom in atoms_obj:
        '%d\t%s\t%f' % (atom.index, atom.symbol, atom.magmom)
    