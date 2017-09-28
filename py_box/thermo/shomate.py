# -*- coding: utf-8 -*-
"""
Created on Thu Dec 01 10:51:14 2016

@author: Jonathan Lym
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import py_box.constants as c

class shomate(object):
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
        t = T/1000        
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

    def get_GoRT(self, T, verbose = True):
        """Calculates the dimensionless free energy (i.e. G/RT) given a temperature."""
        HoRT = self.get_HoRT(T, verbose)
        SoR = self.get_SoR(T, verbose)
        return HoRT-SoR

class shomates(object):
    """
    An object that stores a list of shomate objects.
    """    
    def __init__(self, shomates = [], verbose = True):
        self._shomates = []
        self.verbose = verbose
        for shomate in shomates:
            self.append(shomate)
            
    def append(self, shomate):
        self._shomates.append(shomate)

    def extend(self, shomates):
        self._shomates.extend(shomates)

    def index(self, symbol):
        for i, shomate in enumerate(self):
            if shomate.symbol == symbol:
                return i

    def remove(self, symbol):
        for i, shomate in enumerate(self):
            if shomate.symbol == symbol:
                self._shomates.pop(i)

    def __len__(self):
        return len(self._shomates)
    
    def __setitem__(self, index, shomate):
        self._shomates[index] = shomate

    def __getitem__(self, index):
        return self._shomates[index]

    def _assign_verbose(self, verbose):
        """
        Function that allows verbose input to override the object's verbose for the duration of the calling function.
        """
        if verbose is None:
            return self.verbose
        else:
            return verbose
        
    def print_symbols(self):
        """
        Prints a summary of the thermdat list.
        """
        print "Index\tSymbol"
        for i, shomate in enumerate(self):
            print "%d\t%s" % (i, shomate.symbol)

    def print_rxn(self, stoich_vector):
        """Converts a stoichiometric vector into a stoichiometric reaction."""
        rxn = ''
        react = np.where(np.array(stoich_vector) < 0)[0]
        prod = np.where(np.array(stoich_vector) > 0)[0]
        j = 0    
        #Dealing with reactants
        for i in react:
            rxn += self._add_species_to_rxn(stoich_vector[i], self[i].symbol)
            j += 1        
            if j < len(react):
                rxn += ' + '
        rxn += ' --> '
        j = 0
        #Dealing with products
        for i in prod:
            rxn += self._add_species_to_rxn(stoich_vector[i], self[i].symbol)
            j += 1
            if j < len(react):
                rxn += ' + '
        return rxn

    def _add_species_to_rxn(self, coeff, symbol):
        if abs(coeff) == 1:
            rxn = '%s' % symbol
        else:
            rxn = '%d%s' % (abs(coeff), symbol)
        return rxn
        
    def rxn_to_stoich(self, reaction):
        stoich_vector = [0]*len(self)
        arrow_options = ['-->', '=', '<>', '<->' '<-->', '->', '==']
        reaction_groups = reaction.split(' ') 
        #Find the arrow to determine what are products and reactants
        for arrow_option in arrow_options:
            try:
                arrow_i = reaction_groups.index(arrow_option)
            except ValueError:
                continue
            else:
                break
        else:
            print "Warning. Arrow not found!"
            
        #Go through the groups to determine the species and their stoichiometric coefficients
        for i, reaction_group in enumerate(reaction_groups):
            if i != arrow_i and reaction_group != '+':
                #Determine if reactant or product
                if i < arrow_i:
                    side = -1 #Reactant
                else:
                    side = 1 #Product
                
                #Separate the stoichiometric factor (if any) from the species
                coeff_str = ''             
                for j, character in enumerate(reaction_group):
                    if character.isdigit():
                        coeff_str += character
                    else:
                        coeff_num = 1
                        species = reaction_group[j:]
                        break
                else:
                    print "Warning. Reaction group %s made of only numbers." % reaction_group
                if coeff_str != '':
                    coeff_num = int(float(coeff_str))
                
                #Find the symbol to update the stoichiometric vector
                for j, shomate_species in enumerate(self):
                    if shomate_species.symbol == species:
                        stoich_vector[j] += coeff_num*side
                        break
                else:
                    print "Warning. Could not find the species %s." % species
        return stoich_vector



    def get_HoRT_rxn(self, T = c.T0, reaction = None, stoich_vector = None, H_correction = False, verbose = False):
        """Finds the dimensionless enthalpy of reaction."""
        try:
            iter(T)
        except TypeError:
            HoRT_rxn = 0.
            list_calc = False
        else:
            HoRT_rxn = [0.]*len(T)
            list_calc = True
    
        if reaction is not None:
            stoich_vector = self.rxn_to_stoich(reaction = reaction)
        elif stoich_vector is not None:
            reaction = self.print_rxn(stoich_vector = stoich_vector)
    
        if verbose and not list_calc:
            print "Reaction: %s" % reaction
            print "T = %.2f K" % T 
            print "-"*(10+len(reaction))
            print "Species\tv\tH/RT"
            
        for i, shomate in zip(stoich_vector, self):
            if i != 0:
                HoRT_species = shomate.get_HoRT(T)
                if verbose and not list_calc:
                    print "%s\t%.2f\t%.2f" % (shomate.symbol, i, HoRT_species)
                HoRT_rxn += HoRT_species*i
        if verbose and not list_calc:
            print "-"*(10+len(reaction))
            print "Total H/RT\t\t%.2f" % HoRT_rxn
            print "-"*(10+len(reaction))
        return HoRT_rxn
    
       
    def get_SoR_rxn(self, T = c.T0, reaction = None, stoich_vector = None, verbose = False):
        """Finds the dimensionless entropy of reaction."""
        try:
            iter(T)
        except TypeError:
            SoR_rxn = 0.
            list_calc = False
        else:
            SoR_rxn = [0.]*len(T)    
            list_calc = True
    
        if reaction is not None:
            stoich_vector = self.rxn_to_stoich(reaction = reaction)
        elif stoich_vector is not None:
            reaction = self.print_rxn(stoich_vector = stoich_vector)
    
        if verbose and not list_calc:
            print "Reaction: %s" % reaction
            print "T = %.2f K" % T 
            print "-"*(10+len(reaction))
            print "Species\tv\tS/R"
    
        for i, shomate in zip(stoich_vector, self):
            if i != 0:
                SoR_species = shomate.get_SoR(T)
                if verbose and not list_calc:
                    print "%s\t%.2f\t%.2f" % (shomate.symbol, i, SoR_species)
                SoR_rxn += SoR_species*i
        if verbose and not list_calc:
            print "-"*(10+len(reaction))
            print "Total S/R\t\t%.2f" % SoR_rxn
            print "-"*(10+len(reaction))
        return SoR_rxn
    
    def get_GoRT_rxn(self, T, reaction = None, stoich_vector = None, verbose = True):
        """Finds the dimensionless free energy of reaction."""
        try:
            iter(T)
        except TypeError:
            GoRT_rxn = 0.
            list_calc = False
        else:
            GoRT_rxn = [0.]*len(T)
            list_calc = True
    
        if reaction is not None:
            stoich_vector = self.rxn_to_stoich(reaction = reaction)
        elif stoich_vector is not None:
            reaction = self.print_rxn(stoich_vector = stoich_vector)
    
        if verbose and not list_calc:
            HoRT_rxn = 0.
            SoR_rxn = 0.
    
            print "Reaction: %s" % reaction
            print "T = %.2f K" % T 
            print "-"*(10+len(reaction))
            print "Species\tv\tH/RT\tS/R\tG/RT"
            
        for i, shomate in zip(stoich_vector, self):
            if i != 0:
                HoRT_species = shomate.get_HoRT(T)
                SoR_species = shomate.get_SoR(T)
                GoRT_species = HoRT_species - SoR_species
                if verbose and not list_calc:
                    print "%s\t%.2f\t%.2f\t%.2f\t%.2f" % (shomate.symbol, i, HoRT_species, SoR_species, GoRT_species)
                    HoRT_rxn += HoRT_species*i
                    SoR_rxn += SoR_species*i
                GoRT_rxn += GoRT_species*i
        if verbose and not list_calc:
            print "-"*(10+len(reaction))
            print "Total\t\t%.2f\t%.2f\t%.2f" % (HoRT_rxn, SoR_rxn, GoRT_rxn)
            print "-"*(10+len(reaction))
        return GoRT_rxn
    
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
    
    
def read_shomate(shomate_path, verbose = True):
    shomates_obj = shomates()
    with open(shomate_path, 'r') as shomate_file:
        for line in shomate_file:
            #If the line is not a comment
            if line[0] != '!':
                data = line.split(',')
                if verbose:
                    print "Importing %s" % data[0]
                shomates_obj.append(shomate(symbol = data[0], T_low = float(data[1]), T_high = float(data[2]), a = np.array([float(i) for i in data[3:] if i != '' and i != '\n'])))
    return shomates_obj
