# -*- coding: utf-8 -*-
"""
Created on Thu Dec 01 10:51:14 2016

@author: Jonathan Lym
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import py_box.constants as c

class Shomate(object):
    """
    Contains the Shomate polynomials and corresponding temperature ranges for a species.
    symbol - string to identify the species the object belongs to.
    T_low - Low temperature (K)
    T_high - High temperature (K)
    a - [8x1] Numpy array that holds the Shomate coefficients used between T_low and T_high.
    """
    def __init__(self, symbol, T_low = 0, T_high = 0, a = None):
        self.symbol = symbol
        self.T_low = T_low
        self.T_high = T_high
        if a is None:
            a = np.array(8*[0.])
        self.a = a

    def __str__(self):
        return str(self.__dict__)
    
    def __eq__(self, other):
        if all([all(self.a == other.a), 
                self.T_high == other.T_high,
                self.T_low == other.T_low]): 
            return True
        else:
            return False
        
    def _get_single_CpoR(self, T, verbose = True):
        """Calculates the heat capacity at constant pressure (i.e. Cp/R) given a temperature."""
        t = T/1000.        
        T_arr = np.array([1, t, t ** 2, t ** 3, (1/t) ** 2, 0, 0, 0])
        if verbose:
            if T < self.T_low:
                print "Warning. Input temperature (%f) lower than T_low (%f) for species %s" % (T, self.T_low, self.symbol)
            elif T > self.T_high:
                print "Warning. Input temperature (%f) higher than T_high (%f) for species %s" % (T, self.T_high, self.symbol)
        return np.dot(T_arr, self.a)/c.R('J/mol/K')
    
    def get_CpoR(self, T, verbose = True):
        """Calculates the heat capacity at constant pressure (i.e. Cp/R) given a temperature or a list of temperatures."""
        try: 
            T_val = iter(T)
        except TypeError:
            #Single value T
            CpoR = self._get_single_CpoR(T, verbose)
        else:
            #List value T
            CpoR = np.array([0.]*len(T))
            for i, T_val in enumerate(T):
                CpoR[i] = self._get_single_CpoR(T_val, verbose)            
        return CpoR

    def _get_single_HoRT(self, T, H_correction = False, verbose = True):
        """Calculates the dimensionless enthalpy (i.e. H/RT) given a temperature."""
        t = T/1000.
        if H_correction:
            T_arr = np.array([t, t ** 2 / 2, t ** 3 / 3, t ** 4 / 4, -1/t, 1., 0., 1.])
        else:            
            T_arr = np.array([t, t ** 2 / 2, t ** 3 / 3, t ** 4 / 4, -1/t, 1., 0., 0.])
        if verbose:
            if T < self.T_low:
                print "Warning. Input temperature (%f) lower than T_low (%f) for species %s" % (T, self.T_low, self.symbol)
            elif T > self.T_high:
                print "Warning. Input temperature (%f) higher than T_high (%f) for species %s" % (T, self.T_high, self.symbol)
        return np.dot(T_arr, self.a)/(c.R('kJ/mol/K')*T)

    def get_HoRT(self, T, H_correction = False, verbose = True):
        """Calculates the dimensionless enthalpy at constant pressure (i.e. H/RT) given a temperature or a list of temperatures."""
        try: 
            T_val = iter(T)
        except TypeError:
            #Single value T
            HoRT = self._get_single_HoRT(T, H_correction = H_correction, verbose = verbose)
        else:
            #List value T
            HoRT = np.array([0.]*len(T))
            for i, T_val in enumerate(T):
                HoRT[i] = self._get_single_HoRT(T_val, H_correction = H_correction, verbose = verbose)
        return HoRT

    def _get_single_SoR(self, T, verbose = True):
        """Calculates the dimensionless entropy (i.e. S/R) given a temperature."""
        t = T/1000.        
        T_arr = np.array([np.log(t), t, t ** 2 / 2., t ** 3 / 3., -0.5 * ( 1 / t ) ** 2, 0., 1., 0.])
        if verbose:
            if T < self.T_low:
                print "Warning. Input temperature (%f) lower than T_low (%f) for species %s" % (T, self.T_low, self.symbol)
            elif T > self.T_high:
                print "Warning. Input temperature (%f) higher than T_high (%f) for species %s" % (T, self.T_high, self.symbol)
        return np.dot(T_arr, self.a)/c.R('J/mol/K')

    def get_SoR(self, T, verbose = True):
        """Calculates the dimensionless entropy at constant pressure (i.e. S/R) given a temperature or a list of temperatures."""
        try: 
            T_val = iter(T)
        except TypeError:
            #Single value T
            SoR = self._get_single_SoR(T, verbose)
        else:
            #List value T
            SoR = np.array([0.]*len(T))
            for i, T_val in enumerate(T):
                SoR[i] = self._get_single_SoR(T_val, verbose)   
        return SoR

    def get_GoRT(self, T, H_correction = False, verbose = True):
        """Calculates the dimensionless free energy (i.e. G/RT) given a temperature."""
        HoRT = self.get_HoRT(T, H_correction = H_correction, verbose = verbose)
        SoR = self.get_SoR(T, verbose)
        return HoRT-SoR
    
def _shomate_Cp(t, A, B, C, D, E):
    return A + B*t + C * t**2 + D * t ** 3 + E / t**2
    
def read_fund_csv(symbol, csv_path, print_graph = False):
    """Reads a csv file for data that will be fed into function fit_shomate_species."""
    H0_S0_read = False
    T = []
    Cp = []
    
    print "Reading from file: %s" % csv_path
    with open(csv_path, 'r') as csv_file:
        for line in csv_file:
            if line[0] != '!':
                #Split data and remove unnecessary characters
                data = line.split(',')
                T.append(float(data[0]))
                Cp.append(float(data[1]))
                if len(data) > 2 and not H0_S0_read:
                    H0 = float(data[2])
                    S0 = float(data[3])
                    H0_S0_read = True                

    shomate = fit_shomate_species(symbol, T, Cp, H0, S0)
    if print_graph:
        T_range = np.linspace(shomate.T_low, shomate.T_high)
        Cp_fit = shomate.get_CpoR(T_range)*c.R('J/mol/K')
        H_fit = shomate.get_HoRT(T_range)*T_range*c.R('kJ/mol/K')
        S_fit = shomate.get_SoR(T_range)*c.R('J/mol/K')
        plt.figure()
        plt.subplot(311)
        plt.plot(T, Cp, 'ro', T_range, Cp_fit, 'b-')
        plt.legend(['NIST Data', 'Fit'])
        plt.xlabel('Temperature (K)')
        plt.ylabel('Cp (J/mol/K)')
        
        plt.subplot(312)
        plt.plot(T_range, H_fit)
        plt.xlabel('Temperature (K)')
        plt.ylabel('H (kJ/mol/K)')

        plt.subplot(313)
        plt.plot(T_range, S_fit)
        plt.xlabel('Temperature (K)')
        plt.ylabel('S (J/mol/K)')
    return shomate
    
def fit_shomate_species(symbol, T, Cp, H0, S0):
    """Derives the shomate species from fitting the heat capacity (J/mol/K) and temperature (K) data and including the formation of enthalpy (kJ/mol) and entropy (J/mol/K)."""
    T_low = min(T)
    T_high = max(T)
    t = np.array(T)/1000.
    [a, pcov] = curve_fit(_shomate_Cp, t, np.array(Cp))
    a = list(a)
    a.extend([0., 0., 0.])
    shomate_instance = shomate(symbol = symbol, T_low = T_low, T_high = T_high, a = np.array(a))
    a6 = H0 - shomate_instance.get_HoRT(c.T0)*c.R('kJ/mol/K')*c.T0
    a7 = S0 - shomate_instance.get_SoR(c.T0)*c.R('J/mol/K')
    shomate_instance.a[5] = a6
    shomate_instance.a[6] = a7
    return shomate_instance