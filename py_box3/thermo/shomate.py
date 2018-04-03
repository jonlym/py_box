# -*- coding: utf-8 -*-
"""
Created on Thu Dec 01 10:51:14 2016

@author: Jonathan Lym
"""
import numpy as np
from scipy.optimize import curve_fit
from py_box3 import interpolate
import py_box3.constants as c
import json

class Shomate(object):
    """
    Contains the Shomate polynomials and corresponding temperature ranges for a species.
    symbol - string to identify the species the object belongs to.
    T_low - Low temperature (K)
    T_high - High temperature (K)
    a - [8x1] Numpy array that holds the Shomate coefficients used between T_low and T_high.
    """
    def __init__(self, symbol, T_low = 0, T_high = 0, a = None, phases = None, elements = None):
        self.symbol = symbol
        if a is None:
            a = np.array(8*[0.])
        if type(phases) is str:
            phase_name = phase
        else:
            phase_name = symbol
        if phases is None:
            phases = {symbol: Phase(symbol = phase_name, T_low = T_low, T_high = T_high, a = a)}
        self.phases = phases      
        self.elements = elements


    def _get_single_CpoR(self, T, verbose = True, phase = None):
        if phase is None or phase is 'stable':
            phase = self.get_stable_phase(T = T)
        elif phase is 'temperature':
            phase = self.get_T_phase(T = T)
        return self.phases[phase].get_CpoR(T = T, verbose = verbose)
    
    def get_CpoR(self, T, verbose = True, phase = None):
        """Calculates the heat capacity at constant pressure (i.e. Cp/R) given a temperature or a list of temperatures."""
        try: 
            T_val = iter(T)
        except TypeError:
            #Single value T
            CpoR = self._get_single_CpoR(T = T, verbose = verbose, phase = phase)
        else:
            #List value T
            CpoR = np.array([0.]*len(T))
            for i, T_val in enumerate(T):
                CpoR[i] = self._get_single_CpoR(T = T_val, verbose = verbose, phase = phase)
        return CpoR

    def _get_single_HoRT(self, T, H_correction = False, verbose = True, phase = None):
        if phase is None or phase is 'stable':
            phase = self.get_stable_phase(T = T)
        elif phase is 'temperature':
            phase = self.get_T_phase(T = T)
        return self.phases[phase].get_HoRT(T = T, verbose = verbose, H_correction = H_correction)

    def get_HoRT(self, T, H_correction = False, verbose = True, phase = None):
        """Calculates the dimensionless enthalpy at constant pressure (i.e. H/RT) given a temperature or a list of temperatures."""
        try: 
            T_val = iter(T)
        except:
            #Single value T
            HoRT = self._get_single_HoRT(T = T, H_correction = H_correction, verbose = verbose, phase = phase)
        else:
            #List value T
            HoRT = np.array([0.]*len(T))
            for i, T_val in enumerate(T):
                HoRT[i] = self._get_single_HoRT(T = T_val, H_correction = H_correction, verbose = verbose, phase = phase)
        return HoRT

    def _get_single_SoR(self, T, verbose = True, phase = None):
        """Calculates the dimensionless entropy (i.e. S/R) given a temperature."""
        if phase is None or phase is 'stable':
            phase = self.get_stable_phase(T = T)
        elif phase is 'temperature':
            phase = self.get_T_phase(T = T)
        return self.phases[phase].get_SoR(T = T, verbose = verbose)

    def get_SoR(self, T, verbose = True, phase = None):
        """Calculates the dimensionless entropy at constant pressure (i.e. S/R) given a temperature or a list of temperatures."""
        try: 
            T_val = iter(T)
        except:
            #Single value T
            SoR = self._get_single_SoR(T = T, verbose = verbose, phase = phase)
        else:
            #List value T
            SoR = np.array([0.]*len(T))
            for i, T_val in enumerate(T):
                SoR[i] = self._get_single_SoR(T = T_val, verbose = verbose, phase = phase)
        return SoR

    def get_GoRT(self, T, H_correction = False, verbose = True, phase = None):
        """Calculates the dimensionless free energy (i.e. G/RT) given a temperature."""
        HoRT = self.get_HoRT(T = T, H_correction = H_correction, verbose = verbose, phase = phase)
        SoR = self.get_SoR(T = T, verbose = verbose, phase = phase)
        return HoRT-SoR

    def get_stable_phase(self, T):
        min_phase = None
        min_GoRT = float("inf")
        for name, phase in self.phases.items():
            GoRT = phase.get_GoRT(T = T, verbose = False)
            if GoRT < min_GoRT:
                min_phase = name
                min_GoRT = GoRT
        return min_phase

    def get_T_phase(self, T):
        min_phase = None
        min_GoRT = float("inf")
        for name, phase in self.phases.items():
            GoRT = phase.get_GoRT(T = T, verbose = False)            
            if T <= phase.T_high and T >= phase.T_low and GoRT < min_GoRT:
                min_phase = name
                min_GoRT = GoRT
        return min_phase

    def plot_thermo(self, T_low, T_high, units = None, phase = None):
        """
        Plots the heat capacity, enthalpy and entropy in the temperature range specified.
        The units for the plots can be specified by using R
        """
        import matplotlib.pyplot as plt

        T = np.linspace(T_low, T_high)
        Cp = self.get_CpoR(T, phase = phase)
        H = self.get_HoRT(T, phase = phase)
        S = self.get_SoR(T, phase = phase)
        if units is not None:
            Cp = Cp * c.R(units)
            H = H * c.R(units) * T
            S = S * c.R(units)

        plt.figure()
        plt.subplot(311)
        plt.plot(T, Cp, 'r-')
        if units is None:
            plt.ylabel('Cp/R')
        else:
            plt.ylabel('Cp (%s)' % units)
        plt.xlabel('T (K)')
        plt.title('Plots for %s using shomate polynomials.' % self.symbol)

        plt.subplot(312)
        plt.plot(T, H, 'b-')
        if units is None:
            plt.ylabel('H/RT')
        else:
            plt.ylabel('H (%s)' % units.replace('/K', ''))
        plt.xlabel('T (K)')

        #Entropy graph
        plt.subplot(313)
        plt.plot(T, S, 'k-')
        if units is None:
            plt.ylabel('S/R')
        else:
            plt.ylabel('S (%s)' % units)
        plt.xlabel('T (K)')

    def save_to_json(self, filename = None):
        if filename is None:
            filename = '{}.txt'.format(self.symbol)
        with open(filename, 'w') as f_ptr:
            json.dump(self, f_ptr)

    @classmethod
    def fit_shomate_species(cls, symbol, T, Cp, H0, S0, T_ref = c.T0('K'), elements = None):
        """Derives the shomate species from fitting the heat capacity (J/mol/K) and temperature (K) data and including the formation of enthalpy (kJ/mol) and entropy (J/mol/K)."""
        T_low = min(T)
        T_high = max(T)
        t = np.array(T)/1000.
        [a, pcov] = curve_fit(_shomate_Cp, t, np.array(Cp))
        a = np.append(a, [0., 0., 0.])
        a[5] = H0 - _get_HoRT(T = T_ref, a = a)*c.R('kJ/mol/K')*T_ref
        a[6] = S0 - _get_SoR(T = T_ref, a = a)*c.R('J/mol/K')
        a[7] = - _get_HoRT(T = c.T0('K'), a = a)*c.R('kJ/mol/K')*c.T0('K')
        return cls(symbol = symbol, T_low = T_low, T_high = T_high, a = np.array(a), elements = elements)

    @staticmethod
    def read_fund_csv(symbol, csv_path, print_graph = False):
        """Reads a csv file for data that will be fed into function fit_shomate_species."""
        H0_S0_read = False
        T = []
        Cp = []
        
        print(("Reading from file: %s" % csv_path))
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

        shomate = Shomate.fit_shomate_species(symbol, T, Cp, H0, S0)
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
   
def _shomate_Cp(t, A, B, C, D, E):
    return A + B*t + C * t**2 + D * t ** 3 + E / t**2 

def generate_Cp_data(Ts, Cps, n = 100):
    Ts_out = np.linspace(min(Ts), max(Ts), n)
    Cps_out = np.zeros(shape = n)

    for i, T in enumerate(Ts_out):
        #Finds nearby values to use for interpolation
        low_i = np.argwhere(T >= Ts)[-1]
        high_i = np.argwhere(T <= Ts)[0]
        low_Cp = Cps[low_i]
        high_Cp = Cps[high_i]
        low_T = Ts[low_i]
        high_T = Ts[high_i]

        Cps_out[i] = interpolate(x_low = low_T, x_high = high_T, y_low = low_Cp, y_high = high_Cp, x = T)
    return (Ts_out, Cps_out)

def _get_HoRT(T, a, H_correction = False):
    t = T/1000.
    T_arr = np.array([t, t ** 2 / 2, t ** 3 / 3, t ** 4 / 4, -1/t, 1., 0., 1.])
    return np.dot(T_arr, a)/(c.R('kJ/mol/K')*T)

def _get_SoR(T, a):
    t = T/1000.        
    T_arr = np.array([np.log(t), t, t ** 2 / 2., t ** 3 / 3., -0.5 * ( 1 / t ) ** 2, 0., 1., 0.])
    return np.dot(T_arr, a)/c.R('J/mol/K')



class Phase(object):
    def __init__(self, symbol = '', T_low = 0., T_high = 0., a = None):
        self.symbol = symbol
        self.T_low = T_low
        self.T_high = T_high
        self.a = a

    def get_CpoR(self, T, verbose = True):
        """Calculates the heat capacity at constant pressure (i.e. Cp/R) given a temperature."""
        t = T/1000.        
        T_arr = np.array([1, t, t ** 2, t ** 3, (1/t) ** 2, 0, 0, 0])
        if verbose:
            if T < self.T_low:
                print(("Warning. Input temperature (%f) lower than T_low (%f) for species %s" % (T, self.T_low, self.symbol)))
            elif T > self.T_high:
                print(("Warning. Input temperature (%f) higher than T_high (%f) for species %s" % (T, self.T_high, self.symbol)))
        return np.dot(T_arr, self.a)/c.R('J/mol/K')
    
    def get_HoRT(self, T, H_correction = False, verbose = True):
        """Calculates the dimensionless enthalpy (i.e. H/RT) given a temperature."""
        t = T/1000.
        if H_correction:
            T_arr = np.array([t, t ** 2 / 2, t ** 3 / 3, t ** 4 / 4, -1/t, 1., 0., 1.])
        else:            
            T_arr = np.array([t, t ** 2 / 2, t ** 3 / 3, t ** 4 / 4, -1/t, 1., 0., 0.])
        if verbose:
            if T < self.T_low:
                print(("Warning. Input temperature (%f) lower than T_low (%f) for species %s" % (T, self.T_low, self.symbol)))
            elif T > self.T_high:
                print(("Warning. Input temperature (%f) higher than T_high (%f) for species %s" % (T, self.T_high, self.symbol)))
        return np.dot(T_arr, self.a)/(c.R('kJ/mol/K')*T)

    def get_SoR(self, T, verbose = True):
        """Calculates the dimensionless entropy (i.e. S/R) given a temperature."""
        t = T/1000.        
        T_arr = np.array([np.log(t), t, t ** 2 / 2., t ** 3 / 3., -0.5 * ( 1 / t ) ** 2, 0., 1., 0.])
        if verbose:
            if T < self.T_low:
                print(("Warning. Input temperature (%f) lower than T_low (%f) for species %s" % (T, self.T_low, self.symbol)))
            elif T > self.T_high:
                print(("Warning. Input temperature (%f) higher than T_high (%f) for species %s" % (T, self.T_high, self.symbol)))
        return np.dot(T_arr, self.a)/c.R('J/mol/K')

    def get_GoRT(self, T, H_correction = False, verbose = True):
        """Calculates the dimensionless free energy (i.e. G/RT) given a temperature."""
        HoRT = self.get_HoRT(T, H_correction = H_correction, verbose = verbose)
        SoR = self.get_SoR(T, verbose)
        return HoRT-SoR