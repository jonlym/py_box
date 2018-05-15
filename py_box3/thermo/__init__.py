# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 14:42:40 2016

@author: Jonathan Lym
"""

import os
from ase.io import read
import numpy as np
import matplotlib.pyplot as plt
from pprint import pprint
from py_box3.thermo.thermdat import Thermdat
from py_box3.thermo.thermdats import Thermdats
from py_box3.thermo.nasa import Nasa
import py_box3.constants as c

def dft_to_thermdat(
    input_path, 
    ref_path = 'thermdat_ref.csv',
    T_low = 300.,
    T_high = 800.,
    write_files = True,
    out_path = 'thermdat',
    verbose = False,
    warn = False,
    freq_cut_off = 0.,
    pressure = 1.,
    add_gas_species = True,
    gas_path = 'thermdat_gas'):

    """
    Convert DFT energies and frequencies to NASA polynomials that can be written as a thermdat file
    Attributes
    ----------
        input_path - string
            Path to .csv file that contains the DFT information of all species to be written in thermdat
        ref_path - string
            Path to .csv file that contains the gas-phase reference information
        T_low - float
            Lower temperature bound for the NASA polynomials
        T_high - float
            Higher temperature bound for the NASA polynomials
        write_files - boolean
            Whether or not the thermdat file should be written
        out_path - string
            Path where the thermdat file will be written
        verbose - boolean
            Whether or not more detailed output should be given
        warn - boolean
            Whether or not warnings should be given
        freq_cut_off - float
            Cut off frequency. If vibrational frequencies (in 1/cm) are below this value, they are ignored
            since they are assumed to be frustrated rotational or translational modes
        pressure - float
            Pressure (in atmospheres) to generate the NASA polynomials. Usually this value is left at 1 atm
        add_gas_species - boolean
            Whether or not gas-phase species from gas_path should be written to the thermdat file
        gas_path - string
            Path to thermdat file that contains the gas-phase species to add
    Returns
    -------
        thermdat - Thermdat object
            The successfully converted thermdat object
    """
    T_range = np.linspace(T_low, T_high, T_high-T_low)
    n_T = len(T_range)    

    print("Opening reference file: %s" % ref_path)
    (fund_energies, fund_CHON, rm_list) = read_ref(ref_path, verbose = verbose, warn = warn)

    #Read the species to be processed
    print("Opening input file: %s" % input_path)
    thermdats_dft = read_freq(input_path, verbose = verbose, freq_cut_off = freq_cut_off, warn = warn)

    for i, thermdat_dft in enumerate(thermdats_dft):
        print("-"*10)
        print("Processing %s..." % thermdat_dft.symbol)
        thermdat_dft.nasa = Nasa(symbol = thermdat_dft.symbol, T_low = min(T_range), T_high = max(T_range))    

        print("Calculating heat capacities.")
        CpoR = thermdat_dft.get_CpoR(T_range)
        
        print("Calculating DFT enthalpy of formation")
        if thermdat_dft.is_gas:
            H0 = thermdat_dft.IdealGasThermo.get_enthalpy(temperature = c.T0('K'), verbose = verbose)
        else:
            if np.sum(thermdat_dft.HarmonicThermo.vib_energies) == 0:
                H0 = thermdat_dft.HarmonicThermo.potentialenergy
            else:
                H0 = thermdat_dft.HarmonicThermo.get_internal_energy(temperature = c.T0('K'), verbose = verbose)
        print("Adjusting enthalpy to reference")
        CHON = np.delete(thermdat_dft.CHON, rm_list)
        TransformCHON = np.dot(CHON,np.linalg.inv(fund_CHON))
        HoRT0 = (H0 - np.dot(TransformCHON, fund_energies))/(c.T0('K')*c.kb('eV/K'))

        print("Calculating DFT enthalpy of formation")
        if thermdat_dft.is_gas:
            S0 = thermdat_dft.IdealGasThermo.get_entropy(temperature = c.T0('K'), pressure = pressure*c.convert_unit(from_ = 'atm', to = 'Pa'), verbose = verbose)
        else:
            if np.sum(thermdat_dft.HarmonicThermo.vib_energies) == 0:
                S0 = 0.
            else:
                S0 = thermdat_dft.HarmonicThermo.get_entropy(temperature = c.T0('K'), verbose = verbose)
        SoR0 = S0/c.kb('eV/K')
        
        print("Calculating NASA polynomials")
        thermdat_dft.nasa.fit_NASA(T_range, CpoR, HoRT0, SoR0)

    #Add gas-phase species
    if add_gas_species:
        print("Opening gas phase file: {}".format(gas_path))
        thermdats_gas = Thermdats.from_thermdat(thermdat_path = gas_path, verbose = verbose, warn = warn)
        print("Attaching gas species")
        thermdats_dft.extend(thermdats_gas)

    if write_files:
        print("Writing to thermdat")
        thermdats_dft.write_thermdat(out_path, verbose = False)

    return thermdats_dft

def read_freq(path, freq_cut_off = 0, verbose = True, warn = True):
    """
    Reads the csv file that has the frequencies of the species.
    The CSV file must be of the format:
    Symbol, is_gas, #C, #H, #O, #N, H0, potentialenergy, geometry, symmetry, frequencies

    Parameters
    ----------
        path - string
            Path to the .csv file with the input information
        freq_cut_off - float
            Frequency cut off. If frequencies inputted are less than this value, they are ignored
            since they are assumed to be frustrated rotational modes
        verbose - boolean
            Whether or not more detailed information should be printed
        warn - boolean
            Whether or not warnings should be printed

    Returns
    -------
        thermdat - Thermdat object
            Thermdat object without NASA polynomials
    """
    if verbose:
        print(("Reading from file: %s" % path))
    thermdats = Thermdats()
    with open(path, 'r') as freq_file:
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
                    spin = float(data[10])
                #atoms
                if data[11] == '':
                    atoms = None
                else:
                    atoms = read(data[11])
                #vib_freq
                vib_freq = np.array([float(i) for i in data[12:] if i != '' and i != '\n'])
                if verbose:
                    print(("Ignoring frequencies below %f cm^-1" % freq_cut_off))
                vib_freq = vib_freq[vib_freq >= freq_cut_off]
                if verbose:
                    print(("Importing %s" % symbol))
                thermdats.append(Thermdat(symbol = symbol,
                                              is_gas = is_gas,
                                              CHON = CHON,
                                              H0 = H0,
                                              potentialenergy = potentialenergy,
                                              geometry = geometry,
                                              symmetrynumber = symmetrynumber,
                                              spin = spin,
                                              atoms = atoms,                                            
                                              vib_freq = vib_freq,
                                              verbose = verbose,
                                              warn = warn))
    return thermdats                               

def read_ref(path, freq_cut_off = 0., verbose = True, warn = True):
    """
    Read the reference file to adjust DFT energies to real energies

    Parameters
    ----------
        path - string
            Path to the .csv file that will be read
        freq_cut_off - float
            Frequency cut off. If frequencies inputted are less than this value, they are ignored
            since they are assumed to be frustrated rotational modes
        verbose - boolean
            Whether or not more detailed information should be printed
        warn - boolean
            Whether or not warnings should be printed
    Returns
    -------
        energies_fund - (N,) ndarray
            Fundamental energies of the reference species. Used as a correction factor
        fund_species_CHON - (N, 4) ndarray
            Number of carbon, hydrogen, oxygen and nitrogen in the fundamental species.
        rm_list - list
            Indices to be removed from the fund_species_CHON since the element is not
            present in any of the reference species
    """

    #Read the reference files
    species_ref = read_freq(path, freq_cut_off = freq_cut_off, verbose = verbose, warn = warn)
    
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
        H0_dft[i] = species.IdealGasThermo.get_enthalpy(c.T0('K'), verbose = False)
        H0_exp[i] = species.H0*c.convert_unit(from_ ='kJ/mol', to = 'eV/molecule') #species.H0 from literature

        #Sets up the CHON matrix
        species_CHON[i, :] = species.CHON
        e_sum_list += species.CHON #Used to determine if elements are not needed. e.g. O, N not needed for hydrocarbons    

    #Cleans up the CHON matrix if necessary    
    for i, (element, e_sum) in reversed(list(enumerate(zip(e_list, e_sum_list)))):
        if e_sum == 0:
            if verbose:
                print(('Element %s not needed' % element))
            species_CHON = np.delete(species_CHON, i, 1)
            fund_species_CHON = np.delete(fund_species_CHON, i, 0)
            fund_species_CHON = np.delete(fund_species_CHON, i, 1)
            rm_list.append(i)
    n_e = len(fund_species_CHON)
    if n_e != n_s and verbose:
        print(("Warning: Mismatch between number of elements (%d) and number of reference species." % (n_e, n_s)))
                    
    #Finds the fundamental energy
    H0_fund = np.array([0.]*n_e) #Formation enthalpy for fundamental species. Usually 0
    ref_species_from_fund = np.dot(species_CHON, np.linalg.inv(fund_species_CHON)) #Matrix to go from fundamental species --> reference species
    H0_rxn = H0_exp - np.dot(ref_species_from_fund, H0_fund) #Since H0_fund = 0, this is H0_exp
    energies_fund = np.linalg.solve(ref_species_from_fund, (H0_dft-H0_rxn)) #Energy of fundamental species. Used as correction factor
    return (energies_fund, fund_species_CHON, rm_list)