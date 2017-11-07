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
import matplotlib.pyplot as plt
import py_box.constants as c
import warnings
from py_box.thermo.nasa import Nasa


class thermdat(object):
    """
    Class that holds data related to writing thermdat.
    Required parameters:
    symbol - String that will be written on the first line for 
             the species.
    is_gas - Boolean that is True for gas species and False for surface species. 
             This value will be written as G for True or S for False in the 
             thermdat file.
             Gas species will have an ase.thermochemistry.IdealGasThermdat 
             object while surface species will have an ase.thermochemistry.
             HarmonicThermo object. Each will have additional parameters

    Optional parameters:
    vib_freq - Numpy array that holds the vibrational frequencies in cm^-1.
               If you've been given frequencies in Hz, divide by speed of 
               light in cm/s (29979245800) before passing.
    CHON - List that holds the number of C, H, O and N in that order
    aux_element - String that contains the element name if not C, N, O or H.
    nasa - NASA object that holds temperature ranges and NASA polynomials.
    geometry - String only required for gas species. Options include:
                monatomic                   
                linear
                nonlinear
    potentialenergy - Potential energy obtained from DFT in eV
    symmetrynumber - Symmetry number required for gas species. Input can either
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
    atoms - ASE atoms object only required for gas species
    spin - Integer for the total electronic spin. 
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
        
class thermdats(object):
    """
    An object that stores a list of thermdat objects.
    """    
    def __init__(self, thermdats = [], verbose = True):
        self._thermdats = list(thermdats)
        self.verbose = verbose
            
    def append(self, thermdat):
        self._thermdats.append(thermdat)

    def extend(self, thermdats):
        self._thermdats.extend(thermdats)

    def index(self, symbol):
        for i, thermdat in enumerate(self._thermdats):
            if thermdat.symbol == symbol:
                return i

    def remove(self, symbol):
        for i, thermdat in enumerate(self._thermdats):
            if thermdat.symbol == symbol:
                self._thermdats.pop(i)

    def __len__(self):
        return len(self._thermdats)
    
    def __setitem__(self, index, thermdat):
        self._thermdats[index] = thermdat

    def __getitem__(self, index):
        return self._thermdats[index]
        
    def __str__(self):
        symbols = ''
        for thermdat in self:
            symbols += '%s, ' % thermdat.symbol
        symbols = symbols[:-2]
        return symbols

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
        for i, thermdat in enumerate(self):
            print "%d\t%s" % (i, thermdat.symbol)

    def write_thermdat(self, thermdat_path, verbose = None):
        """
        Writes the thermdat file using the list of thermdat species
        """
        verbose = self._assign_verbose(verbose)
        with open(thermdat_path, 'w') as thermdat_file:
            if verbose:
                print "Writing thermdat to path: %s" % thermdat_path
            thermdat_file.write('THERMO ALL\n       100       500      1500\n')
            
            float_string = '%.8E'    
            for thermdat in self:
                if verbose:
                    print 'Writing species: %s' % thermdat.symbol
                self._write_line1(thermdat_file, thermdat)
                self._write_line2(thermdat_file, thermdat, float_string)
                self._write_line3(thermdat_file, thermdat, float_string)
                self._write_line4(thermdat_file, thermdat, float_string)
            thermdat_file.write('END')    
            
    def _write_line1(self, thermdat_file, species):
        """
        Writes the first line of the thermdat file
        """
        line1_pos = [
            24, # C
            28, # C#
            29, # H
            33, # H #
            34, # O
            38, # O #
            39, # N
            43, # N #
            44, # Phase
            45, # T_low
            55, # T_high
            65, # T_mid
            79] # Line num
    
        if species.is_gas:
            phase = 'G'
        else:
            phase = 'S'
    
        #Adjusts the position based on the number of CHON        
        for i in range(len(species.CHON)):
            line1_pos[i*2+1] = line1_pos[i*2+1] - len(str(int(species.CHON[i]))) + 1
    
        line1_fields = [
            species.symbol,
            'C',
            '%d' % species.CHON[0],
            'H',
            '%d' % species.CHON[1],
            'O',
            '%d' % species.CHON[2],
            'N',
            '%d' % species.CHON[3],
            phase,
            '%.1f' % species.nasa.T_low,
            '%.1f' % species.nasa.T_high,
            '%.1f' % species.nasa.T_mid,
        ]
        line = ''
        for pos, field in zip(line1_pos, line1_fields):
            line += field
            line = self._insert_space(pos, line)
        line += '1\n'
        thermdat_file.write(line)
    
    def _write_line2(self, thermdat_file, species, float_string):
        """
        Writes the second line of the thermdat file
        """
        line = ''
        for i in range(5):
            a = species.nasa.a_high[i]
            if a >= 0:
                line += ' '            
            line += float_string % a
        line += '    2\n'
        thermdat_file.write(line)
    
    def _write_line3(self, thermdat_file, species, float_string):
        """
        Writes the third line of the thermdat file
        """
        line = ''
        for i in range(5):
            if i < 2:
                a = species.nasa.a_high[i+5]
            else:
                a = species.nasa.a_low[i-2]
            if a >= 0:
                line += ' '
            line += float_string % a
        line += '    3\n'
        thermdat_file.write(line)
    
    def _write_line4(self, thermdat_file, species, float_string):
        """
        Writes the fourth line of the thermdat file
        """
        line = ''
        for i in range(3,7):
            a = species.nasa.a_low[i]
            if a >= 0:
                line += ' '
            line += float_string % a
        line += '                   4\n'
        thermdat_file.write(line)
    
    def _insert_space(self, end_index, string):
        """
        Inserts the number of spaces required given the string and the position of
        the next non-blank field.
        """
        string += ' ' * (end_index - len(string))
        return string

    def get_HoRT_rxn(self, T = c.T0, reaction = None, stoich_vector = None, verbose = None):
        """Finds the dimensionless enthalpy of reaction."""
        verbose = self._assign_verbose(verbose)
        try:
            iter(T)
        except TypeError:
            HoRT_rxn = 0.
            list_calc = False
        else:
            HoRT_rxn = [0.]*len(T)
            list_calc = True
    
        if reaction is not None:
            stoich_vector = self.reaction_to_stoich(reaction = reaction)
        elif stoich_vector is not None:
            reaction = self.print_reaction(stoich_vector = stoich_vector)
    
        if verbose and not list_calc:
            print "Reaction: %s" % reaction
            print "T = %.2f K" % T 
            print "-"*(10+len(reaction))
            print "Species\tv\tH/RT"
            
        for i, thermdat in zip(stoich_vector, self):
            if i != 0:
                HoRT_species = thermdat.nasa.get_HoRT(T)
                if verbose and not list_calc:
                    print "%s\t%.2f\t%.2f" % (thermdat.symbol, i, HoRT_species)
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
            stoich_vector = self.reaction_to_stoich(reaction = reaction)
        elif stoich_vector is not None:
            reaction = self.print_reaction(stoich_vector = stoich_vector)
    
        if verbose and not list_calc:
            print "Reaction: %s" % reaction
            print "T = %.2f K" % T 
            print "-"*(10+len(reaction))
            print "Species\tv\tS/R"
    
        for i, thermdat_species in zip(stoich_vector, self):
            if i != 0:
                SoR_species = thermdat_species.nasa.get_SoR(T)
                if verbose and not list_calc:
                    print "%s\t%.2f\t%.2f" % (thermdat_species.symbol, i, SoR_species)
                SoR_rxn += SoR_species*i
        if verbose and not list_calc:
            print "-"*(10+len(reaction))
            print "Total S/R\t\t%.2f" % SoR_rxn
            print "-"*(10+len(reaction))
        return SoR_rxn
    
    def get_GoRT_rxn(self, T, reaction = None, stoich_vector = None, verbose = False):
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
            stoich_vector = self.reaction_to_stoich(reaction = reaction)
        elif stoich_vector is not None:
            reaction = self.print_reaction(stoich_vector = stoich_vector)
    
        if verbose and not list_calc:
            HoRT_rxn = 0.
            SoR_rxn = 0.
    
            print "Reaction: %s" % reaction
            print "T = %.2f K" % T 
            print "-"*(10+len(reaction))
            print "Species\tv\tH/RT\tS/R\tG/RT"
            
        for i, thermdat in zip(stoich_vector, self):
            if i != 0:
                HoRT_species = thermdat.nasa.get_HoRT(T)
                SoR_species = thermdat.nasa.get_SoR(T)
                GoRT_species = HoRT_species - SoR_species
                if verbose and not list_calc:
                    print "%s\t%.2f\t%.2f\t%.2f\t%.2f" % (thermdat.symbol, i, HoRT_species, SoR_species, GoRT_species)
                    HoRT_rxn += HoRT_species*i
                    SoR_rxn += SoR_species*i
                GoRT_rxn += GoRT_species*i
        if verbose and not list_calc:
            print "-"*(10+len(reaction))
            print "Total\t\t%.2f\t%.2f\t%.2f" % (HoRT_rxn, SoR_rxn, GoRT_rxn)
            print "-"*(10+len(reaction))
        return GoRT_rxn
        
    def write_to_excel(self, file_name = 'thermdat.xlsx'):
        """Takes a thermdat list and writes it to file_name."""
        print "Creating %s workbook" % file_name
        book = xlwt.Workbook()
        #sh = book.get_sheet('Sheet 1')
        sh = book.add_sheet('Sheet 1')
    
        print "Adding row headings"
        row_headings = [' ', 'T low', 'T_mid', 'T_high', 'Low Range', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', ' ', 'High Range', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7']
        col = 0
        for row, row_heading in enumerate(row_headings):
            sh.write(row, col, row_heading)
    
        for col, thermdat_species in enumerate(self):
            print "Processing %s" % thermdat_species.symbol
            col_offset = 1
            row = 0
            sh.write(row, col+col_offset, thermdat_species.symbol)
            sh.write(row+1, col+col_offset, thermdat_species.nasa.T_low)
            sh.write(row+2, col+col_offset, thermdat_species.nasa.T_mid)
            sh.write(row+3, col+col_offset, thermdat_species.nasa.T_high)
            
            for row, (a_low, a_high) in enumerate(zip(thermdat_species.nasa.a_low, thermdat_species.nasa.a_high)):
                a_low_offset = 5
                a_high_offset = 14
                sh.write(row + a_low_offset, col+col_offset, a_low)
                sh.write(row + a_high_offset, col+col_offset, a_high)    
        book.save(file_name)
        
    def print_reaction(self, stoich_vector):
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
        rxn += ' = '
        j = 0
        #Dealing with products
        for i in prod:
            rxn += self._add_species_to_rxn(stoich_vector[i], self[i].symbol)
            j += 1
            if j < len(prod):
                rxn += ' + '
        return rxn
    
    def _add_species_to_rxn(self, coeff, symbol):
        if abs(coeff) == 1:
            rxn = '%s' % symbol
        else:
            rxn = '%d%s' % (abs(coeff), symbol)
        return rxn
        
    def reaction_to_stoich(self, reaction):
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
            warnings.warn("Arrow not found!")
            
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
                    warnings.warn("Reaction group %s made of only numbers." % reaction_group)
                if coeff_str != '':
                    coeff_num = int(float(coeff_str))
                
                #Find the symbol to update the stoichiometric vector
                for j, thermdat_species in enumerate(self):
                    if thermdat_species.symbol == species:
                        stoich_vector[j] += coeff_num*side
                        break
                else:
                    warnings.warn("Could not find the species %s." % species)
        return stoich_vector


def read_thermdat(thermdat_path, verbose = True, warn = True):
    """
    Reads the thermdat file and returns an list of molecule objects with the
    symbols, temperature ranges and CHON fields.
    """
    if verbose:
        print "Reading thermdat file: %s" % thermdat_path
    thermdats_obj = thermdats()
    #Possible atoms 
    atoms = ['C', 'H', 'O', 'N']
    site_type = 1
    #Positions at which coefficients start in lines 2, 3, 4
    pos_index = [0, 15, 30, 45, 60]
    with open(thermdat_path, 'r') as thermdat_file:
        for line in thermdat_file:
            if len(line) > 2:
                #Finds the line number
                line_num = line[-3:]
                if '1' in line_num:
                    #Encountered new species
                    symbol_index = line.find(' ')
                    symbol = line[0:symbol_index]
                    if verbose:
                        print "Importing %s" % symbol
                    buf= line[symbol_index:-1]
                    CHON = []
                    for atom in atoms:
                        CHON.append(_get_CHON_value(buf, symbol, atom, verbose, warn))
                    try:
                        phase = re.search("(G|S)", buf).group(0)
                    except AttributeError:
                        warnings.warn("Phase not found. Assuming surface phase.")
                        is_gas = False
                    else:
                        if 'G' in phase:
                            is_gas = True
                        elif 'S' in phase:
                            is_gas = False
                    T_range = _get_T_values(line, symbol)
                elif '2' in line_num:
                    a_high = []
                    for i in pos_index:
                        a_high.append(float(line[i:i+15]))
                elif '3' in line_num:
                    a_low  = []
                    for i in pos_index:
                        if len(a_high) < 7:
                            a_high.append(float(line[i:i+15]))
                        else:
                            a_low.append(float(line[i:i+15]))
                elif '4' in line_num:
                    for i in pos_index:
                        if len(a_low) < 7:
                            a_low.append(float(line[i:(i+15)]))
                    nasa = Nasa(T_low = T_range[0],
                                T_mid = T_range[2], 
                                T_high = T_range[1], 
                                a_low = np.array(a_low), 
                                a_high = np.array(a_high),
                                verbose = verbose,
                                symbol = symbol)
                    thermdats_obj.append(thermdat(symbol = symbol, 
                                                 CHON = CHON, 
                                                 is_gas = is_gas, 
                                                 nasa = nasa,
                                                 site_type = site_type,
                                                 verbose = verbose,
                                                  warn = warn))
    return thermdats_obj

def _get_CHON_value(string, symbol, atom, verbose = True, warn = True):
    """Returns the number of C, H, O or N in an atom based on the string inputted."""
    #Looks for a pattern which starts with the atom type (e.g. C), has spaces
    #and then ends with a digit    
    pattern = '%s +\d+' % atom  
    try:
        buf = re.search(pattern, string).group(0)
    except AttributeError:
        if warn:
            warnings.warn("Unable to find %s atom in species %s. Returning 0." % (atom, symbol))
        return 0
    buf = re.search('\d+', buf).group(0)
    return int(float(buf))

def _get_T_values(string, symbol, verbose = True):
    """Returns T_low, T_mid and T_high given a string inputted."""
    #Looks for three floating numbers separated by a space
    pattern = "[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)? +[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)? +[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?"
    try:
        buf = re.search(pattern, string).group(0)
    except AttributeError:
        if verbose:
            print "Warning: Unable to find temperature limits in species %s. Returning 0s." % (symbol)
        return [0, 0, 0]
    T_limits = re.split(" +", buf)
    T_out = []
    for T in T_limits:
        T_out.append(float(T))
    if len(T_out) < 3:
        if verbose:
            print "Warning: Not all temperatures found for species %s. Returning 0s for unfound values." % symbol
        while len(T_out) < 3:
            T_out.append(0)
    return T_out
