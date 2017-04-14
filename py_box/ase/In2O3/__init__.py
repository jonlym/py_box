# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 17:12:09 2017

@author: Jonathan Lym
"""

import heapq

index_dict = {'In1': [26, 27],
              'In2': [24, 25],
              'In3': [15, 14],
              'In4': [22, 23],
              'O1': [75, 74],
              'O2': [76, 77],
              'O3': [46, 47],
              'O4': [73, 72],
              'O5': [43, 42],
              'O6': [68, 69]}

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
