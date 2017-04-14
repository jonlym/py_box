# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 16:38:29 2016

@author: Jonathan Lym
"""

import numpy as np
from datetime import datetime, timedelta
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
        
