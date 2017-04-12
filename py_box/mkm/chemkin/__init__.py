# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 14:57:39 2016

@author: Jonathan Lym
"""

from py_box.constants import T0, convert_unit
from ase.io import read
from py_box.mkm.chemkin.thermdat import thermdat, thermdats
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

def read_ref(path, verbose = True):
    #Read the reference files
    species_ref = read_freq(path, verbose = verbose)
    
    #Prepare the fundamental species
    e_list = ['C', 'H', 'O', 'N']
    rm_list = []
    fund_species_CHON = np.array([[1, 0, 0, 0],
                                 [0, 2, 0, 0],
                                 [0, 0, 2, 0],
                                 [0, 0, 0, 2]])
    
    n_s = len(species_ref)
    H0_dft = np.array([0.]*n_s)
    H0_exp = np.array([0.]*n_s)
    species_CHON = np.array([[0]*4]*n_s)
    e_sum_list = np.array([0.]*4)
    
    for i, species in enumerate(species_ref):
        #Calculate properties at T0
        H0_dft[i] = species.IdealGasThermo.get_enthalpy(T0, verbose = False)
        H0_exp[i] = species.H0*convert_unit(from_ ='kJ/mol', to = 'eV/molecule') #species.H0 from literature

        #Sets up the CHON matrix
        species_CHON[i, :] = species.CHON
        e_sum_list += species.CHON #Used to determine if elements are not needed. e.g. O, N not needed for hydrocarbons    

    #Cleans up the CHON matrix if necessary    
    for i, (element, e_sum) in reversed(list(enumerate(zip(e_list, e_sum_list)))):
        if e_sum == 0:
            if verbose:
                print 'Element %s not needed' % element
            species_CHON = np.delete(species_CHON, i, 1)
            fund_species_CHON = np.delete(fund_species_CHON, i, 0)
            fund_species_CHON = np.delete(fund_species_CHON, i, 1)
            rm_list.append(i)
    n_e = len(fund_species_CHON)
    if n_e != n_s and verbose:
        print "Warning: Mismatch between number of elements (%d) and number of reference species." % (n_e, n_s)
                    
    #Finds the fundamental energy
    H0_fund = np.array([0.]*n_e) #Formation enthalpy for fundamental species. Usually 0
    ref_species_from_fund = np.dot(species_CHON, np.linalg.inv(fund_species_CHON)) #Matrix to go from fundamental species --> reference species
    H0_rxn = H0_exp - np.dot(ref_species_from_fund, H0_fund) #Since H0_fund = 0, this is H0_exp
    energies_fund = np.linalg.solve(ref_species_from_fund, (H0_dft-H0_rxn)) #Energy of fundamental species. Used as correction factor
    return (energies_fund, fund_species_CHON, rm_list)