# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 16:11:06 2017

@author: Jonathan Lym
"""

from py_box import any_alpha
from ase.io import read
from ase.io.bader import attach_charges
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

class bader:
    """Stores the bader charge"""
    def __init__(self, charges = [], min_distances = [], atomic_volumes = [], atoms = None):
        self.charges = charges
        self.min_distances = min_distances
        self.atomic_volumes = atomic_volumes
        self.atoms = atoms

    def append(self, atom, charge = np.nan, min_distance = np.nan, atomic_volume = np.nan,):
        """Appends data. Note that an atom object must be supplied to ensure consistency."""
        self.charges.append(charge)
        self.min_distancs.append(min_distance)
        self.atomic_volumes.append(atomic_volume)
        self.atoms.append(atom)

    @classmethod
    def from_ACF_dat(cls, ACF_path = 'ACF.dat', atoms = None):
        """Reads the bader charge from the ACF.dat file"""
        charges = []
        min_distances = []
        atomic_volumes = []
        with open(ACF_path, 'r') as ACF_file:
            for line in ACF_file:
                #If there are no alphabetical characters and does not contain a line
                if (not any_alpha(line) and '--' not in line):
                    data = [np.float(x) for x in line.split(' ') if x != '']
                    charges.append(data[4])
                    min_distances.append(data[5])
                    atomic_volumes.append(data[6])
        return cls(charges = charges, min_distances = min_distances, atomic_volumes = atomic_volumes, atoms = atoms)

    def __str__(self):
        lines = []
        lines.append("Element[index]  Charges  Min Distance  Atomic Volume\n")
        lines.append("----------------------------------------------------\n")
        for atom, charge, min_distance, atomic_volume in zip(self.atoms, self.charges, self.min_distances, self.atomic_volumes):
            atom_info = '{}[{}]'.format(atom.symbol, atom.index)
            pad_length = len('Element[index]  ') - len(atom_info)
            lines.append('{}{}{:7.2f}  {:12.2f}{:13.2f}\n'.format(atom_info, ' '*pad_length, charge, min_distance, atomic_volume))
        return ''.join(lines)

def print_bader(dir_path = '.', atoms_file = 'CONTCAR'):
    #Prints the bader charges using the atoms_file (usually CONTCAR) and ACF.dat
    atoms = read('%s/%s' % (dir_path, atoms_file))
    attach_charges(atoms, '%s/ACF.dat' % dir_path)
    buf = []    
    print "i\tSymbol\tCharge\tx\ty\tz"
    print "-"*30
    for atom in atoms:
        print "%d\t%s\t%f\t%f\t%f\t%f" % (atom.index, atom.symbol, atom.charge, atom.x, atom.y, atom.z)
        buf = buf.append((atom.index, atom.symbol, atom.charge, atom.x, atom.y, atom.z))        
    dtype = [('i', int), ('symbol', str) ('x', float), ('y', float), ('z', float)]
    data = np.array(buf, dtype = dtype)
    return data

        
def compare_bader(dir1, dir2):
    #Uses the CONTCARs located in dir1 and dir2 and takes the difference between the two charges
    atoms1 = read('%s/CONTCAR' % dir1)
    attach_charges(atoms1, '%s/ACF.dat' % dir1)
    atoms2 = read('%s/CONTCAR' % dir2)
    attach_charges(atoms2, '%s/ACF.dat' % dir2)
    max_len = max([len(atoms1, atoms2)])

    len1 = len(atom1)
    len2 = len(atom2)

    if len1 >= len2:
        print "Using properties of atom1."
        max_len = len1
        atoms_long = atoms1
        atoms_short = atoms2
    else:
        print "Using properties of atom2 since it is longer."
        max_len = len2
        atoms_long = atoms2
        atoms_short = atoms1
    
    print "i\tSymbol\tCharge 1\t Charge 2\t ΔCharge"
    for i in range(max_len):
        try:
            atom_short = atoms_short[i]
        except IndexError: 
            print "%d\t%s\t%f\t%f\t%f" % (atom1.index, atom1.symbol, atom1.charge, atom2.charge, atom2.charge - atom1.charge)
        else:
            print "%d\t%s\t%f\t%f\t%f" % (atom1.index, atom1.symbol, atom1.charge, atom2.charge, atom2.charge - atom1.charge)            
        buf = buf.append((atom1.index, atom1.symbol, atom1.charge, atom1.x, atom1.y, atom1.z))        
        if atom1.symbol != atom2.symbol:
            print "Symbol mismatch!"
    dtype = [('i', int), ('symbol', str) ('x', float), ('y', float), ('z', float)]
    data = np.array(buf, dtype = dtype)    
    return data
