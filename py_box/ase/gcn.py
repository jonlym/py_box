# -*- coding: utf-8 -*-
"""
Created on Sun May 07 17:00:00 2017

@author: Jon Lym
"""

import os, warnings
import numpy as np
from ase.geometry import find_mic
from py_box import get_unique_list
try:
    import pandas as pd
except:
    pass
try:
    import xlsxwriter
except:
    pass

class GCN(object):
    """
    Contains parameters relating to the coordination number of the atoms.
    Parameters:
        atoms - ASE atoms object.
        bulk - ASE atoms object of the bulk material.
        bulk_CN - Dict
            Bulk coordination number for each element.
        atom_radii - Dict
            Radii of each atom. The function import_atom_radii may be used to import from the reference file.
        scale - Float
            Used to increase (>1) or decrease (<1) the effective atomic radii.
        exceptions - Dict
            Holds the list of elements that cannot be coordinated.
            e.g. exceptions = {'Cu': ['H', 'O']}
            Copper cannot be coordinated to hydrogen or oxygen.
        source_dict - Dict
            The key are the elements present in the atoms object and the value can be:
                'Empirical'
                'Calculated'
                'van der Waals'
                'Covalent (single bond)'
                'Covalent (triple bond)'
                'Metallic'
                Custom values by entering a float
            e.g. source_dict = {'Cu': 'Metallic',
                                'H': 'Empirical'}
    """
    
    def __init__(self, atoms, bulk, scale = 1., exceptions = {}, source_dict = None):
        self.atoms = atoms
        self.bulk = bulk
        self.bulk_CN = {}
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
            self.get_bulk_CN()
            self.calc_CNs()

        #Goes through neighbors and determines their coordination number
        for i, (atom, CN) in enumerate(zip(self.atoms, self.CNs)):
            try:
                self.GCNs[i] = CN/self.bulk_CN[atom.symbol]
            except KeyError:
                continue
        # for i in range(len(self.atoms)):
        #     if self.CNs[i] == 0:
        #         self.GCNs[i] = 0
        #     else:
        #         neighbors_CNs = np.zeros(len(self.neighbors[i]))
        #         for j, k in enumerate(self.neighbors[i]):
        #             neighbors_CNs[j] = self.CNs[k]
        #         self.GCNs[i] = np.sum(neighbors_CNs)/np.max(neighbors_CNs)

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

    def get_bulk_CN(self):
        """
        Calculates the coordination number and the list of indices that the
        bulk atoms are coordinated to.
        """

        #Duplicate the bulk atoms object to account for periodic images
        offset = 13 * len(self.bulk)
        bulk_ext = self.bulk.copy() * (3, 3, 3)
        for i in range(offset, (offset + len(self.bulk))):
            symbol = bulk_ext[i].symbol

            CN = 0
            for j in range(len(bulk_ext)):
                #If not the same index, radii overlap (weighted by scale) and the elements are not in the element list
                if (i != j and
                    bulk_ext.get_distance(i, j) <= ((self.atom_radii[bulk_ext[i].symbol] + self.atom_radii[bulk_ext[j].symbol]) * self.scale) and
                    bulk_ext[j].symbol not in self.exceptions[bulk_ext[i].symbol]):
                    #Add a neighbor
                    CN += 1

            #Add record to dictionary
            if bulk_ext[i].symbol not in self.bulk_CN:
                self.bulk_CN[symbol] = 0
            if CN > self.bulk_CN[symbol]:
                self.bulk_CN[symbol] = CN

#Radius of atoms. Source: https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
atom_radii_dict = {
    'H': {'Empirical':              0.25,
          'Calculated':             0.53,
          'van der Waals':          1.20,
          'Covalent (single bond)': 0.38,
          'Covalent (triple bond)': None,
          'Metallic':               None,
          },
    'He': {'Empirical':              None,
           'Calculated':             0.31,
           'van der Waals':          1.40,
           'Covalent (single bond)': 0.32,
           'Covalent (triple bond)': None,
           'Metallic':               None,
           },
    'Li': {'Empirical':              1.45,
           'Calculated':             1.67,
           'van der Waals':          1.82,
           'Covalent (single bond)': 1.34,
           'Covalent (triple bond)': None,
           'Metallic':               1.52,
           },
    'Be': {'Empirical':              1.05,
           'Calculated':             1.12,
           'van der Waals':          1.53,
           'Covalent (single bond)': 0.9,
           'Covalent (triple bond)': 0.85,
           'Metallic':               1.12,
           },
    'B': {'Empirical':              0.85,
          'Calculated':             0.87,
          'van der Waals':          1.92,
          'Covalent (single bond)': 0.82,
          'Covalent (triple bond)': 0.73,
          'Metallic':               None,
          },
    'C': {'Empirical':              0.70,
          'Calculated':             0.67,
          'van der Waals':          1.70,
          'Covalent (single bond)': 0.77,
          'Covalent (triple bond)': 0.60,
          'Metallic':               None,
          },
    'N': {'Empirical':              0.65,
          'Calculated':             0.56,
          'van der Waals':          1.55,
          'Covalent (single bond)': 0.75,
          'Covalent (triple bond)': 0.54,
          'Metallic':               None,
          },
    'O': {'Empirical':              0.60,
          'Calculated':             0.48,
          'van der Waals':          1.52,
          'Covalent (single bond)': 0.73,
          'Covalent (triple bond)': 0.53,
          'Metallic':               None,
          },
    'F': {'Empirical':              0.50,
          'Calculated':             0.42,
          'van der Waals':          1.47,
          'Covalent (single bond)': 0.71,
          'Covalent (triple bond)': 0.53,
          'Metallic':               None,
          },
    'Ne': {'Empirical':              None,
           'Calculated':             0.38,
           'van der Waals':          1.54,
           'Covalent (single bond)': 0.69,
           'Covalent (triple bond)': None,
           'Metallic':               None,
           },
    'Na': {'Empirical':              1.50,
           'Calculated':             1.45,
           'van der Waals':          1.73,
           'Covalent (single bond)': 1.30,
           'Covalent (triple bond)': 1.27,
           'Metallic':               1.60,
           },
    'Mg': {'Empirical':              1.25,
           'Calculated':             1.18,
           'van der Waals':          1.84,
           'Covalent (single bond)': 1.18,
           'Covalent (triple bond)': 1.11,
           'Metallic':               1.43,
           },
    'Al': {'Empirical':              1.25,
           'Calculated':             1.18,
           'van der Waals':          1.84,
           'Covalent (single bond)': 1.18,
           'Covalent (triple bond)': 1.11,
           'Metallic':               1.43,
           },
    'Si': {'Empirical':              1.10,
           'Calculated':             1.11,
           'van der Waals':          2.10,
           'Covalent (single bond)': 1.11,
           'Covalent (triple bond)': 1.02,
           'Metallic':               None,
           },
    'P': {'Empirical':               1.00,
          'Calculated':              0.98,
          'van der Waals':           1.80,
          'Covalent (single bond)':  1.06,
          'Covalent (triple bond)':  0.94,
          'Metallic':                None,
          },
    'S': {'Empirical':               1.00,
          'Calculated':              0.88,
          'van der Waals':           1.80,
          'Covalent (single bond)':  1.02,
          'Covalent (triple bond)':  0.95,
          'Metallic':                None,
          },
    'Cl': {'Empirical':              1.00,
           'Calculated':             0.79,
           'van der Waals':          1.75,
           'Covalent (single bond)': 0.99,
           'Covalent (triple bond)': 0.93,
           'Metallic':               None,
           },
    'Ar': {'Empirical':              0.71,
           'Calculated':             0.71,
           'van der Waals':          1.88,
           'Covalent (single bond)': 0.97,
           'Covalent (triple bond)': 0.96,
           'Metallic':               None,
           },
    'K': {'Empirical':               2.20,
          'Calculated':              2.43,
          'van der Waals':           2.75,
          'Covalent (single bond)':  1.96,
          'Covalent (triple bond)':  None,
          'Metallic':                2.27,
          },
    'Ca': {'Empirical':              1.80,
           'Calculated':             1.94,
           'van der Waals':          2.31,
           'Covalent (single bond)': 1.74,
           'Covalent (triple bond)': 1.33,
           'Metallic':               1.97,
           },
    'Sc': {'Empirical':              1.60,
           'Calculated':             1.84,
           'van der Waals':          2.11,
           'Covalent (single bond)': 1.44,
           'Covalent (triple bond)': 1.14,
           'Metallic':               1.62,
           },
    'Ti': {'Empirical':              1.40,
           'Calculated':             1.76,
           'van der Waals':          None,
           'Covalent (single bond)': 1.36,
           'Covalent (triple bond)': 1.08,
           'Metallic':               1.47,
           },
    'V': {'Empirical':               1.35,
          'Calculated':              1.71,
          'van der Waals':           None,
          'Covalent (single bond)':  1.25,
          'Covalent (triple bond)':  1.06,
          'Metallic':                1.34,
          },
    'Cr': {'Empirical':              1.40,
           'Calculated':             1.66,
           'van der Waals':          None,
           'Covalent (single bond)': 1.27,
           'Covalent (triple bond)': 1.03,
           'Metallic':               1.28,
           },
    'Mn': {'Empirical':              1.40,
           'Calculated':             1.61,
           'van der Waals':          None,
           'Covalent (single bond)': 1.39,
           'Covalent (triple bond)': 1.03,
           'Metallic':               1.27,
           },
    'Fe': {'Empirical':              1.40,
           'Calculated':             1.56,
           'van der Waals':          None,
           'Covalent (single bond)': 1.25,
           'Covalent (triple bond)': 1.02,
           'Metallic':               1.26,
           },
    'Co': {'Empirical':              1.35,
           'Calculated':             1.52,
           'van der Waals':          None,
           'Covalent (single bond)': 1.26,
           'Covalent (triple bond)': 0.96,
           'Metallic':               1.25,
           },
    'Ni': {'Empirical':              1.35,
           'Calculated':             1.49,
           'van der Waals':          1.63,
           'Covalent (single bond)': 1.21,
           'Covalent (triple bond)': 1.01,
           'Metallic':               1.24,
           },
    'Cu': {'Empirical':              1.35,
           'Calculated':             1.45,
           'van der Waals':          1.40,
           'Covalent (single bond)': 1.38,
           'Covalent (triple bond)': 1.20,
           'Metallic':               1.28,
           },
    'Zn': {'Empirical':              1.35,
           'Calculated':             1.42,
           'van der Waals':          1.39,
           'Covalent (single bond)': 1.31,
           'Covalent (triple bond)': None,
           'Metallic':               1.34,
           },
    'Ga': {'Empirical':              1.30,
           'Calculated':             1.36,
           'van der Waals':          1.87,
           'Covalent (single bond)': 1.26,
           'Covalent (triple bond)': 1.21,
           'Metallic':               1.35,
           },
    'Ge': {'Empirical':              1.25,
           'Calculated':             1.25,
           'van der Waals':          2.11,
           'Covalent (single bond)': 1.22,
           'Covalent (triple bond)': 1.14,
           'Metallic':               None,
           },
    'As': {'Empirical':              1.15,
           'Calculated':             1.14,
           'van der Waals':          1.85,
           'Covalent (single bond)': 1.19,
           'Covalent (triple bond)': 1.06,
           'Metallic':               None,
           },
    'Se': {'Empirical':              1.15,
           'Calculated':             1.03,
           'van der Waals':          1.90,
           'Covalent (single bond)': 1.16,
           'Covalent (triple bond)': 1.07,
           'Metallic':               None,
           },
    'Br': {'Empirical':              1.15,
           'Calculated':             0.94,
           'van der Waals':          1.85,
           'Covalent (single bond)': 1.14,
           'Covalent (triple bond)': 1.10,
           'Metallic':               None,
           },
    'Kr': {'Empirical':              None,
           'Calculated':             0.88,
           'van der Waals':          2.02,
           'Covalent (single bond)': 1.10,
           'Covalent (triple bond)': 1.08,
           'Metallic':               None,
           },
    'Rb': {'Empirical':              2.35,
           'Calculated':             2.65,
           'van der Waals':          3.03,
           'Covalent (single bond)': 2.11,
           'Covalent (triple bond)': None,
           'Metallic':               2.48,
           },
    'Sr': {'Empirical':              2.00,
           'Calculated':             2.19,
           'van der Waals':          2.49,
           'Covalent (single bond)': 1.92,
           'Covalent (triple bond)': 1.39,
           'Metallic':               2.15,
           },
    'Y': {'Empirical':               1.80,
          'Calculated':              2.12,
          'van der Waals':           None,
          'Covalent (single bond)':  1.62,
          'Covalent (triple bond)':  1.24,
          'Metallic':                1.80,
          },
    'Zr': {'Empirical':              1.55,
           'Calculated':             2.06,
           'van der Waals':          None,
           'Covalent (single bond)': 1.48,
           'Covalent (triple bond)': 1.21,
           'Metallic':               1.60,
           },
    'Nb': {'Empirical':              1.45,
           'Calculated':             1.98,
           'van der Waals':          None,
           'Covalent (single bond)': 1.37,
           'Covalent (triple bond)': 1.16,
           'Metallic':               1.46,
           },
    'Mo': {'Empirical':              1.45,
           'Calculated':             1.90,
           'van der Waals':          None,
           'Covalent (single bond)': 1.45,
           'Covalent (triple bond)': 1.13,
           'Metallic':               1.39,
           },
    'Tc': {'Empirical':              1.35,
           'Calculated':             1.83,
           'van der Waals':          None,
           'Covalent (single bond)': 1.56,
           'Covalent (triple bond)': 1.10,
           'Metallic':               1.36,
           },
    'Ru': {'Empirical':              1.30,
           'Calculated':             1.78,
           'van der Waals':          None,
           'Covalent (single bond)': 1.26,
           'Covalent (triple bond)': 1.03,
           'Metallic':               1.34,
           },
    'Rh': {'Empirical':              1.35,
           'Calculated':             1.73,
           'van der Waals':          None,
           'Covalent (single bond)': 1.35,
           'Covalent (triple bond)': 1.06,
           'Metallic':               1.34,
           },
    'Pd': {'Empirical':              1.40,
           'Calculated':             1.69,
           'van der Waals':          1.63,
           'Covalent (single bond)': 1.31,
           'Covalent (triple bond)': 1.12,
           'Metallic':               1.37,
           },
    'Ag': {'Empirical':              1.60,
           'Calculated':             1.65,
           'van der Waals':          1.72,
           'Covalent (single bond)': 1.53,
           'Covalent (triple bond)': 1.37,
           'Metallic':               1.44,
           },
    'Cd': {'Empirical':              1.55,
           'Calculated':             1.61,
           'van der Waals':          1.58,
           'Covalent (single bond)': 1.48,
           'Covalent (triple bond)': None,
           'Metallic':               1.51,
           },
    'In': {'Empirical':              1.55,
           'Calculated':             1.56,
           'van der Waals':          1.93,
           'Covalent (single bond)': 1.44,
           'Covalent (triple bond)': 1.46,
           'Metallic':               1.67,
           },
    'Sn': {'Empirical':              1.45,
           'Calculated':             1.45,
           'van der Waals':          2.17,
           'Covalent (single bond)': 1.41,
           'Covalent (triple bond)': 1.32,
           'Metallic':               None,
           },
    'Sb': {'Empirical':              1.45,
           'Calculated':             1.33,
           'van der Waals':          2.06,
           'Covalent (single bond)': 1.38,
           'Covalent (triple bond)': 1.27,
           'Metallic':               None,
           },
    'Te': {'Empirical':              1.40,
           'Calculated':             1.23,
           'van der Waals':          2.06,
           'Covalent (single bond)': 1.35,
           'Covalent (triple bond)': 1.21,
           'Metallic':               None,
           },
    'I': {'Empirical':               1.40,
          'Calculated':              1.15,
          'van der Waals':           1.98,
          'Covalent (single bond)':  1.33,
          'Covalent (triple bond)':  1.25,
          'Metallic':                None,
          },
    'Xe': {'Empirical':              None,
           'Calculated':             1.08,
           'van der Waals':          2.16,
           'Covalent (single bond)': 1.30,
           'Covalent (triple bond)': 1.22,
           'Metallic':               None,
           },
    'Cs': {'Empirical':              2.60,
           'Calculated':             2.98,
           'van der Waals':          3.43,
           'Covalent (single bond)': 2.25,
           'Covalent (triple bond)': None,
           'Metallic':               2.65,
           },
    'Ba': {'Empirical':              2.15,
           'Calculated':             2.53,
           'van der Waals':          2.68,
           'Covalent (single bond)': 1.98,
           'Covalent (triple bond)': 1.49,
           'Metallic':               2.22,
           },
    'La': {'Empirical':              1.95,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': 1.69,
           'Covalent (triple bond)': 1.39,
           'Metallic':               1.87,
           },
    'Ce': {'Empirical':              1.85,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': 1.31,
           'Metallic':               1.82,
           },
    'Pr': {'Empirical':              1.85,
           'Calculated':             2.47,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': 1.28,
           'Metallic':               1.82,
           },
    'Nd': {'Empirical':              1.85,
           'Calculated':             2.06,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               1.81,
           },
    'Pm': {'Empirical':              1.85,
           'Calculated':             2.05,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               1.83,
           },
    'Sm': {'Empirical':              1.85,
           'Calculated':             2.38,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               1.80,
           },
    'Eu': {'Empirical':              1.85,
           'Calculated':             2.31,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               1.80,
           },
    'Gd': {'Empirical':              1.80,
           'Calculated':             2.33,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': 1.32,
           'Metallic':               1.80,
           },
    'Tb': {'Empirical':              1.75,
           'Calculated':             2.25,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               1.77,
           },
    'Dy': {'Empirical':              1.75,
           'Calculated':             2.28,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               1.78,
           },
    'Ho': {'Empirical':              1.75,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               1.76,
           },
    'Er': {'Empirical':              1.75,
           'Calculated':             2.26,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               1.76,
           },
    'Tm': {'Empirical':              1.75,
           'Calculated':             2.22,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               1.76,
           },
    'Yb': {'Empirical':              1.75,
           'Calculated':             2.22,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               1.76,
           },
    'Lu': {'Empirical':              1.75,
           'Calculated':             2.17,
           'van der Waals':          None,
           'Covalent (single bond)': 1.60,
           'Covalent (triple bond)': 1.31,
           'Metallic':               1.74,
           },
    'Hf': {'Empirical':              1.55,
           'Calculated':             2.08,
           'van der Waals':          None,
           'Covalent (single bond)': 1.50,
           'Covalent (triple bond)': 1.22,
           'Metallic':               1.59,
           },
    'Ta': {'Empirical':              1.45,
           'Calculated':             2.00,
           'van der Waals':          None,
           'Covalent (single bond)': 1.38,
           'Covalent (triple bond)': 1.19,
           'Metallic':               1.46,
           },
    'W': {'Empirical':               1.35,
          'Calculated':              1.93,
          'van der Waals':           None,
          'Covalent (single bond)':  1.46,
          'Covalent (triple bond)':  1.15,
          'Metallic':                1.39,
          },
    'Re': {'Empirical':              1.35,
           'Calculated':             1.88,
           'van der Waals':          None,
           'Covalent (single bond)': 1.59,
           'Covalent (triple bond)': 1.10,
           'Metallic':               1.37,
           },
    'Os': {'Empirical':              1.30,
           'Calculated':             1.85,
           'van der Waals':          None,
           'Covalent (single bond)': 1.28,
           'Covalent (triple bond)': 1.09,
           'Metallic':               1.35,
           },
    'Ir': {'Empirical':              1.35,
           'Calculated':             1.80,
           'van der Waals':          None,
           'Covalent (single bond)': 1.37,
           'Covalent (triple bond)': 1.07,
           'Metallic':               1.36,
           },
    'Pt': {'Empirical':              1.35,
           'Calculated':             1.77,
           'van der Waals':          1.75,
           'Covalent (single bond)': 1.28,
           'Covalent (triple bond)': 1.10,
           'Metallic':               1.39,
           },
    'Au': {'Empirical':              1.35,
           'Calculated':             1.74,
           'van der Waals':          1.66,
           'Covalent (single bond)': 1.44,
           'Covalent (triple bond)': 1.23,
           'Metallic':               1.44,
           },
    'Hg': {'Empirical':              1.50,
           'Calculated':             1.71,
           'van der Waals':          1.55,
           'Covalent (single bond)': 1.49,
           'Covalent (triple bond)': None,
           'Metallic':               1.51,
           },
    'Tl': {'Empirical':              1.90,
           'Calculated':             1.56,
           'van der Waals':          1.96,
           'Covalent (single bond)': 1.48,
           'Covalent (triple bond)': 1.50,
           'Metallic':               1.70,
           },
    'Pb': {'Empirical':              1.80,
           'Calculated':             1.54,
           'van der Waals':          2.02,
           'Covalent (single bond)': 1.47,
           'Covalent (triple bond)': 1.37,
           'Metallic':               None,
           },
    'Bi': {'Empirical':              1.60,
           'Calculated':             1.43,
           'van der Waals':          2.07,
           'Covalent (single bond)': 1.46,
           'Covalent (triple bond)': 1.35,
           'Metallic':               None,
           },
    'Po': {'Empirical':              1.90,
           'Calculated':             1.35,
           'van der Waals':          1.97,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': 1.29,
           'Metallic':               None,
           },
    'At': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          2.02,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': 1.38,
           'Metallic':               None,
           },
    'Rn': {'Empirical':              None,
           'Calculated':             1.20,
           'van der Waals':          2.20,
           'Covalent (single bond)': 1.45,
           'Covalent (triple bond)': 1.33,
           'Metallic':               None,
           },
    'Fr': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          3.48,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               None,
           },
    'Ra': {'Empirical':              2.15,
           'Calculated':             None,
           'van der Waals':          2.83,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': 1.59,
           'Metallic':               None,
           },
    'Ac': {'Empirical':              1.95,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': 1.40,
           'Metallic':               None,
           },
    'Th': {'Empirical':              1.80,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': 1.36,
           'Metallic':               1.79,
           },
    'Pa': {'Empirical':              1.80,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': 1.29,
           'Metallic':               1.63,
           },
    'U': {'Empirical':               1.75,
          'Calculated':              None,
          'van der Waals':           1.86,
          'Covalent (single bond)':  None,
          'Covalent (triple bond)':  1.18,
          'Metallic':                1.56,
          },
    'Np': {'Empirical':              1.75,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': 1.16,
           'Metallic':               1.55,
           },
    'Pu': {'Empirical':              1.75,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               1.59,
           },
    'Am': {'Empirical':              1.75,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               1.73,
           },
    'Cm': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               1.74,
           },
    'Bk': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               1.70,
           },
    'Cf': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               1.86,
           },
    'Es': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               1.86,
           },
    'Fm': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               None,
           },
    'Md': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               None,
           },
    'No': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               None,
           },
    'Lr': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               None,
           },
    'Rf': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': 1.31,
           'Metallic':               None,
           },
    'Db': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': 1.26,
           'Metallic':               None,
           },
    'Sg': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': 1.21,
           'Metallic':               None,
           },
    'Bh': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': 1.19,
           'Metallic':               None,
           },
    'Hs': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': 1.18,
           'Metallic':               None,
           },
    'Mt': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': 1.13,
           'Metallic':               None,
           },
    'Ds': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': 1.12,
           'Metallic':               None,
           },
    'Rg': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': 1.18,
           'Metallic':               None,
           },
    'Cn': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': 1.30,
           'Metallic':               None,
           },
    'Nh': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               None,
           },
    'Fl': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               None,
           },
    'Mc': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               None,
           },
    'Lv': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               None,
           },
    'Ts': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               None,
           },
    'Og': {'Empirical':              None,
           'Calculated':             None,
           'van der Waals':          None,
           'Covalent (single bond)': None,
           'Covalent (triple bond)': None,
           'Metallic':               None,
           }
}
