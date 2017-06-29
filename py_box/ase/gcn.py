# -*- coding: utf-8 -*-
"""
Created on Sun May 07 17:00:00 2017

@author: Jon Lym
"""

import os, warnings
import numpy as np
import pandas as pd
import xlsxwriter

class GCN(object):
    """
    Contains parameters relating to the coordination number of the atoms.
    Parameters:
        atoms - ASE atoms object.
        atom_radii - Dictionary containing the radii of each atom. The function import_atom_radii may be used to import from the reference file.
        scale - A float that is used to increase (>1) or decrease (<1) the effective atomic radii.
        exceptions - Dictionary that holds the list of elements that cannot be coordinated.
                     e.g. exceptions = {'Cu': ['H', 'O']}
                     Copper cannot be coordinated to hydrogen or oxygen.
    """
    
    def __init__(self, atoms, scale = 1., exceptions = {}, source_dict = None):
        self.atoms = atoms
        self.CNs = np.zeros(len(atoms))
        self.neighbors = [[] for i in range(len(self.atoms))]
        self.GCNs = np.zeros(len(atoms))
        self.scale = scale
        self.exceptions = exceptions
        self.atom_radii = {}
        if source_dict is not None:
            self.import_atom_radii(source_dict)

        #If the exception and atom_radii dictionaries are not empty, find the coordination numbers
        if any(self.exceptions) and any(self.atom_radii):
            self.calc_GCNs()

    def import_atom_radii(self, source_dict):
        """
        Assigns the radius in Anstroms to the self.atom_radii. 
        source_dict is a dictionary where the key corresponds to the element and the 
        value can be one of the following:
            'Empirical'
            'Calculated'
            'van der Waals'
            'Covalent (single bond)'
            'Covalent (triple bond)'
            'Metallic'
            Custom values can be specified by entering a float.
        """
        #Import reference file
        current_path = os.getcwd()
        os.chdir(os.environ['PYBOX'])
        atoms_ref = pd.read_excel('atom_radius.xlsx')
        os.chdir(current_path)

        #Assign radius
        for symbol, r_source in source_dict.iteritems():
            #Custom value
            if type(r_source) is float:
                if r_source < 0:
                    warnings.warn('Element %s assigned a negative atomic radius: %f A.' % (symbol, r_source))
                self.atom_radii[symbol] = r_source            
            #Value to be imported
            elif r_source in ['Empirical', 'Calculated', 'van der Waals', 'Covalent (single bond)', 'Covalent (triple bond)', 'Metallic']:
                self.atom_radii[symbol] = atoms_ref[symbol][r_source]
                if self.atom_radii[symbol] < 0.:
                    warnings.warn('Radius type %s for element %s not available in reference table.' % (r_source, symbol))          
            else:
                warnings.warn('Invalid radius type %s for element %s.' % (r_source, symbol))
        
    def calc_CNs(self):
        """
        Calculates the coordination number and the list of indices that the 
        atoms are coordinated to.
        """

        #Duplicate the atoms object to account for periodic images
        offset = 13 * len(self.atoms)
        atoms_ext = self.atoms.copy() * (3, 3, 3)
        for i in range(offset, (offset + len(self.atoms))):
            neighbors = []
            for j in range(len(atoms_ext)):
                #If not the same index, radii overlap (weighted by scale) and the elements are not in the element list
                if (i != j and 
                    atoms_ext.get_distance(i, j) <= ((self.atom_radii[atoms_ext[i].symbol] + self.atom_radii[atoms_ext[j].symbol]) * self.scale) and
                    atoms_ext[j].symbol not in self.exceptions[atoms_ext[i].symbol]):
                    #Add a neighbor
                    self.CNs[i % len(self.atoms)] += 1
                    neighbors.append(j % len(self.atoms))
            self.neighbors[i % len(self.atoms)] = get_unique_list(neighbors)
            
    def calc_GCNs(self, update = True):
        """
        Calculates the generalized coordination number.
        """
        if update:
            self.calc_CNs()

        #Goes through neighbors and determines their coordination number
        for i in range(len(self.atoms)):
            if self.CNs[i] == 0:
                self.GCNs[i] = 0
            else:
                neighbors_CNs = np.zeros(len(self.neighbors[i]))
                for j, k in enumerate(self.neighbors[i]):
                    neighbors_CNs[j] = self.CNs[k]
                self.GCNs[i] = np.sum(neighbors_CNs)/np.max(neighbors_CNs)

    def write_to_excel(self, file_name = 'gcn.xlsx'):
        """
        Writes an excel file containing all the data in GCN.
        """
        workbook = xlsxwriter.Workbook(file_name)
        worksheet = workbook.add_worksheet()
        headers = ['Index', 'Symbol', 'CN', 'GCN', 'Neighbors']
        for j, header in enumerate(headers):
            worksheet.write(0, j, header)

        for i, (atom, CN, GCN, neighbors) in enumerate(zip(self.atoms, self.CNs, self.GCNs, self.neighbors)):
            worksheet.write(i+1, 0, atom.index)
            worksheet.write(i+1, 1, atom.symbol)
            worksheet.write(i+1, 2, CN)
            worksheet.write(i+1, 3, GCN)
            for j, neighbor in enumerate(neighbors):
                worksheet.write(i+1, j+4, '%s%d' % (self.atoms[neighbor].symbol, neighbor))

    def get_generalized_property(self, index, property):
        """
        Given a list corresponding to the indices in the atoms object, calculates the generalized property of the neighbors.
        index: Integer containing the index of the atom.
        property: List containing the properties of the atom (indices aligned with the self.atoms object
        returns: Float of the generalized property
        """
        neighbor_properties = np.zeros(shape = len(self.neighbors[index]))
        for i, neighbor in enumerate(self.neighbors[index]):
            neighbor_properties[i] = property[neighbor]
        return np.sum(neighbor_properties)/np.max(neighbor_properties)

    def get_neighbors_sum_property(self, index, property):
        """
        Given a list corresponding to the indices in the atoms object, calculates the summed property of the neighbors.
        index: Integer containing the index of the atom.
        property: List containing the properties of the atom (indices aligned with the self.atoms object
        returns: Float of the generalized property
        """
        neighbor_properties = np.zeros(shape = len(self.neighbors[index]))
        for i, neighbor in enumerate(self.neighbors[index]):
            neighbor_properties[i] = property[neighbor]
        return np.sum(neighbor_properties)

    def get_neighbors_average_property(self, index, property):
        """
        Given a list corresponding to the indices in the atoms object, calculates the averaged property of the neighbors.
        index: Integer containing the index of the atom.
        property: List containing the properties of the atom (indices aligned with the self.atoms object
        returns: Float of the generalized property
        """
        neighbor_properties = np.zeros(shape = len(self.neighbors[index]))
        for i, neighbor in enumerate(self.neighbors[index]):
            neighbor_properties[i] = property[neighbor]
        return np.average(neighbor_properties)

def get_unique_list(data):
    """
    Given a list, returns a unique list.
    """
    keys = {}
    for item in data:
        keys[item] = 1
    return keys.keys()
