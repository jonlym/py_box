# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 17:12:09 2017

@author: Jonathan Lym
"""

from ase.io import read
import heapq
import numpy as np

index_dict = {'In1': [26, 27],
              'In2': [24, 25],
              'In3': [15, 14],
              'In4': [22, 23],
              'O1':  [75, 74],
              'O2':  [76, 77],
              'O3':  [46, 47],
              'O4':  [73, 72],
              'O5':  [43, 42],
              'O6':  [68, 69],
              'In1A': 26,
              'In1B': 27,
              'In2A': 24,
              'In2B': 25,
              'In3A': 15,
              'In3B': 14,
              'In4A': 22,
              'In4B': 23,
              'O1A':  75,
              'O1B':  74,
              'O2A':  76,
              'O2B':  77,
              'O3A':  46,
              'O3B':  47,
              'O4A':  73,
              'O4B':  72,
              'O5A':  43,
              'O5B':  42,
              'O6A':  68,
              'O6B':  69}

def find_closest_H(atoms, i, n = 1):
    """Finds the 'n' indices of the H atoms closest to atom 'i'."""
    all_H_indices = []
    all_H_dist = []
    for j, atom in enumerate(atoms):
        if atom.symbol == 'H':
            all_H_indices.append(j)
            all_H_dist.append(atoms.get_distance(i, j))
    H_low_dist = heapq.nsmallest(n, all_H_dist)
    H_low_indices = []
    for dist in H_low_dist:
        H_low_indices.append(all_H_indices[all_H_dist.index(dist)])
    return H_low_indices

def get_new_index(atoms_indices, vacancy_indices, start_index = 0):
    """
    Returns the adjusted indices given the vacancies already created on the surface.
    start_index = 1 should be used when the indices are not 0 indexed.
    """
    out_atom_indices = []

    if vacancy_indices is None:
        return atoms_indices
    elif type(vacancy_indices) is int:
        vacancy_indices = [vacancy_indices]

    if type(atoms_indices) is int:
        atom_indices = [atoms_indices]

    for atom_index in atom_indices:
        if atom_index in vacancy_indices:
            out_atom_indices.append(np.nan)
        else:
            offset = start_index - sum([1 for vacancy_index in vacancy_indices if vacancy_index < atom_index])
            out_atom_indices.append(atom_index + offset)

    if len(out_atom_indices) == 1:
        return out_atom_indices[0]
    else:
        return out_atom_indices

def get_In2O3_configuration(n = 0, width = 12):
    indices = [75, 76, 46, 73, 43, 68, 74, 77, 47, 72, 42, 69]
    config = [int(x) for x in ''.split(np.binary_repr(num = n, width = width))]
    del_indices = indices[config]
    atoms = read('/home/work/ccei_biomass/users/jlym/In2O3/unit_cell_relaxed/In2O3_110_clean')
    del atoms[del_indices]
    return atoms
