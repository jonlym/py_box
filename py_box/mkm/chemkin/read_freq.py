# -*- coding: utf-8 -*-
"""
Created on Sun Nov 27 15:58:44 2016

@author: Jon Lym
"""

from ase.io import read
from thermdat import thermdat, thermdats
import numpy as np

def read_freq(freq_file_path, freq_cut_off = 0, verbose = True):
    """
    Reads the csv file that has the frequencies of the species.
    The CSV file must be of the format:
    Symbol, is_gas, #C, #H, #O, #N, H0, potentialenergy, geometry, symmetry, frequencies
    """
    if verbose:
        print "Reading from file: %s" % freq_file_path
    thermdats_obj = thermdats()
    freq_file = open(freq_file_path, 'r')
    with open(freq_file_path, 'r') as freq_file:
        for line in freq_file:
            if line[0] is not '!':
                #Split data and remove unnecessary characters
                data = line.split(',')
                #Sorting data into respective sets
                #Symbol
                symbol = data[0]
                #is_gas            
                is_gas = int(data[1])            
                #CHON            
                CHON = [int(float(i)) for i in data[2:6]]
                if data[6] == '':
                    data[6] = '0'
                H0 = float(data[6])
                #potentialenergy            
                if data[7] == '':
                    data[7] = '0'
                potentialenergy = float(data[7])            
                #geometry
                if data[8] == '':
                    geometry = None
                else:
                    geometry = data[8]
                #symmetrynumber
                if data[9] == '':
                    symmetrynumber = None
                else:
                    symmetrynumber = int(data[9])
                #spin
                if data[10] == '':
                    spin = None
                else:
                    spin = int(data[10])
                #atoms
                if data[11] == '':
                    atoms = None
                else:
                    atoms = read(data[11])
                #vib_freq
                vib_freq = np.array([float(i) for i in data[12:] if i != '' and i != '\n'])
                if verbose:
                    print "Ignoring frequencies below %f cm^-1" % freq_cut_off
                vib_freq = vib_freq[vib_freq >= freq_cut_off]
                if verbose:
                    print "Importing %s" % symbol
                thermdats_obj.append(thermdat(symbol = symbol,
                                              is_gas = is_gas,
                                              CHON = CHON,
                                              H0 = H0,
                                              potentialenergy = potentialenergy,
                                              geometry = geometry,
                                              symmetrynumber = symmetrynumber,
                                              spin = spin,
                                              atoms = atoms,                                            
                                              vib_freq = vib_freq,
                                              verbose = verbose))
    return thermdats_obj                                