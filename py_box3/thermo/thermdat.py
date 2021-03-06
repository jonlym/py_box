# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 14:57:39 2016

@author: Jonathan Lym
"""

from scipy.stats import variation
import re
import numpy as np
import ase.thermochemistry
import xlwt
import py_box3.constants as c
import warnings
from py_box3.thermo.nasa import Nasa


class Thermdat(object):
    """
    Class that holds data related to writing thermdat.

    Required Attributes
    -------------------
        symbol - string
            Will be written on the first line for the species.
        is_gas - boolean 
            True for gas species and False for surface species. 
            This value will be written as G for True or S for False in the 
            thermdat file.
            Gas species will have an ase.thermochemistry.IdealGasThermdat 
            object while surface species will have an ase.thermochemistry.
            HarmonicThermo object. Each will have additional parameters

    Optional Attributes
    -------------------
        vib_freq - (N,) ndarray
            Holds the vibrational frequencies in cm^-1.
            If you've been given frequencies in Hz, divide by speed of 
            light in cm/s (29979245800) before passing.
        CHON - (4,) list
            Holds the number of C, H, O and N in that order
        aux_element - string
            Holds the element name if not C, N, O or H.
        nasa - NASA object 
            Holds temperature ranges and NASA polynomials.
        geometry - string 
            Only required for gas species. Options include:
                monatomic                   
                linear
                nonlinear
        potentialenergy - float
            Potential energy obtained from DFT in eV
        symmetrynumber - int or string
            Symmetry number required for gas species. Input can either
            be an integer or a string of the corresponding point group.
            Supported point groups:
            Point group    symmetry number
            C1             1
            Cs             1
            C2             2
            C2v            2
            C3v            3
            Cinfv          1
            D2h            4
            D3h            6
            D5h            10
            Dinfh          2
            D3d            6
            Td             12
            Oh             24
            See DOI for more details: 10.1007/s00214-007-0328-0
        atoms - ASE atoms object 
            Only required for gas species to calculate rotational modes
        spin - float
            Holds the total electronic spin. 
            0 for molecules in which all electrons are paired
            0.5 for a free radical with a single unpaired electron
            1.0 for a triplet with two unpaired electrons, such as O_2.
    """
    def __init__(self, 
                 symbol, 
                 is_gas, 
                 vib_freq = None, 
                 CHON = None,
                 aux_element = None,
                 nasa = None, 
                 geometry = None,
                 H0 = 0.,
                 potentialenergy = 0.,
                 symmetrynumber = None,
                 atoms = None,
                 spin = None,
                 site_type = 1,
                 verbose = True,
                 warn = True):

        self.verbose = verbose                     
        self.symbol = symbol
        self.is_gas = is_gas
        if not self.is_gas:
            self.site_type = site_type
        if vib_freq is None:
            if warn:
                warnings.warn("%s was not assigned any vibrational frequencies. This is needed for full functionality for thermochemistry." % symbol)
            vib_freq = np.array([0.])
        self._vib_freq = vib_freq
        vib_energies = self.__get_vib_energies()
        self.H0 = H0
        if is_gas:
            if geometry is None:
                geometry = 'monatomic'
                if warn:
                    warnings.warn("%s is a gas species but a geometry was not assigned. \nAssigning the default geometry, %s" % (symbol, geometry))
            if symmetrynumber is None:
                symmetrynumber = 1
                if warn:
                    warnings.warn("%s is a gas species but a symmetry number was not assigned. \nAssigning a default symmetry number, %d" % (symbol, symmetrynumber))
            elif type(symmetrynumber) is str:
                sym_dict = {
                    'C1': 1,
                    'Cs': 1,
                    'C2': 2,
                    'C2v': 2,
                    'C3v': 3,
                    'Cinfv': 1,
                    'D2h': 4,
                    'D3h': 6,
                    'D5h': 10,
                    'Dinfh': 2,
                    'D3d': 6,
                    'Td': 12,
                    'Oh': 24       
                }
                if sym_dict.get(symmetrynumber) is None:
                    symmetrynumber = 1
                    if warn:
                        warnings.warn("Point group not found. Using default symmetry number of %d" % symmetrynumber)
            if atoms is None:
                if warn:
                    warnings.warn("%s is a gas species but an ASE atoms object was not assigned." % symbol)
            if spin is None:
                spin = 0
                if warn:
                    warnings.warn("%s is a gas species but a spin was not assigned. \nAssigning a default value of %d" % (symbol, spin))
            self.IdealGasThermo = ase.thermochemistry.IdealGasThermo(vib_energies = vib_energies, 
                                                                geometry = geometry, 
                                                                potentialenergy = potentialenergy,
                                                                symmetrynumber = symmetrynumber,
                                                                atoms = atoms,
                                                                spin = spin)
        else:
            self.HarmonicThermo = ase.thermochemistry.HarmonicThermo(vib_energies = vib_energies, 
                                                                     potentialenergy = potentialenergy)                                                                
        if CHON is None:
            CHON = 4*[0]
            self.aux_element = aux_element
        self.CHON = CHON
        if nasa is None:
            nasa = Nasa(symbol = symbol)
        self.nasa = nasa

    @property
    def vib_freq(self):
        return self._vib_freq
        
    @vib_freq.setter
    def vib_freq(self, vib_freq):
        self._vib_freq = vib_freq
        if self.is_gas:
            self.IdealGasThermo.vib_energies = self.__get_vib_energies()
        else:
            self.HarmonicThermo.vib_energies = self.__get_vib_energies()            

    def __str__(self):
        return str(self.__dict__)
    
    def __get_vib_energies(self):
        """
        Gets the vibrational energies in eV
        """
        return self._vib_freq * c.c('cm/s') * c.h('eV s')
        
    def _get_single_CpoR(self, T):
        """
        Returns the constant pressure heat capacity.
        """
        if self.is_gas:
            if self.IdealGasThermo.geometry == 'monatomic':
                return 2.5
            elif self.IdealGasThermo.geometry == 'linear':
                CpoR = 2.5
            elif self.IdealGasThermo.geometry == 'nonlinear':
                CpoR = 3.
            else:
                warnings.warn("Gas phase species with invalid geometry. Returning 0.")
                return 0.                
            #Correction between Cv and Cp
            CpoR += 1.
        else:
            CpoR = 0.
        for vib_f in self.vib_freq:               
            vib_T = vib_f*c.c('cm/s')*c.h('J s')/c.kb('J/K')
            CpoR += (vib_T / ( 2. * T ))**2 * 1./(np.sinh(vib_T / ( 2. * T )))**2
        return CpoR

    def get_CpoR(self, T):
        """Calculates the heat capacity at constant pressure (i.e. Cp/R) given a temperature or a list of temperatures."""
        try: 
            T_val = iter(T)
        except TypeError:
            #Single value T
            CpoR = self._get_single_CpoR(T)
        else:
            #List value T
            CpoR = np.array([0.]*len(T))
            for i, T_val in enumerate(T):
                CpoR[i] = self._get_single_CpoR(T_val)            
        return CpoR