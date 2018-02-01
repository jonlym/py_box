import numpy as np
from py_box.thermo.shomate import Shomate
from py_box import constants as c

class Shomates(object):
    """
    An object that stores a list of shomate objects.
    """    
    def __init__(self, shomates = [], verbose = True):
        self._shomates = list(shomates)
        self.verbose = verbose
            
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
        print("Index\tSymbol")
        for i, shomate in enumerate(self):
            print(("%d\t%s" % (i, shomate.symbol)))

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
            print("Warning. Arrow not found!")
            
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
                    print(("Warning. Reaction group %s made of only numbers." % reaction_group))
                if coeff_str != '':
                    coeff_num = int(float(coeff_str))
                
                #Find the symbol to update the stoichiometric vector
                for j, shomate_species in enumerate(self):
                    if shomate_species.symbol == species:
                        stoich_vector[j] += coeff_num*side
                        break
                else:
                    print(("Warning. Could not find the species %s." % species))
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
            print(("Reaction: %s" % reaction))
            print(("T = %.2f K" % T)) 
            print(("-"*(10+len(reaction))))
            print("Species\tv\tH/RT")
            
        for i, shomate in zip(stoich_vector, self):
            if i != 0:
                HoRT_species = shomate.get_HoRT(T)
                if verbose and not list_calc:
                    print(("%s\t%.2f\t%.2f" % (shomate.symbol, i, HoRT_species)))
                HoRT_rxn += HoRT_species*i
        if verbose and not list_calc:
            print(("-"*(10+len(reaction))))
            print(("Total H/RT\t\t%.2f" % HoRT_rxn))
            print(("-"*(10+len(reaction))))
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
            print(("Reaction: %s" % reaction))
            print(("T = %.2f K" % T)) 
            print(("-"*(10+len(reaction))))
            print("Species\tv\tS/R")
    
        for i, shomate in zip(stoich_vector, self):
            if i != 0:
                SoR_species = shomate.get_SoR(T)
                if verbose and not list_calc:
                    print(("%s\t%.2f\t%.2f" % (shomate.symbol, i, SoR_species)))
                SoR_rxn += SoR_species*i
        if verbose and not list_calc:
            print(("-"*(10+len(reaction))))
            print(("Total S/R\t\t%.2f" % SoR_rxn))
            print(("-"*(10+len(reaction))))
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
    
            print(("Reaction: %s" % reaction))
            print(("T = %.2f K" % T)) 
            print(("-"*(10+len(reaction))))
            print("Species\tv\tH/RT\tS/R\tG/RT")
            
        for i, shomate in zip(stoich_vector, self):
            if i != 0:
                HoRT_species = shomate.get_HoRT(T)
                SoR_species = shomate.get_SoR(T)
                GoRT_species = HoRT_species - SoR_species
                if verbose and not list_calc:
                    print(("%s\t%.2f\t%.2f\t%.2f\t%.2f" % (shomate.symbol, i, HoRT_species, SoR_species, GoRT_species)))
                    HoRT_rxn += HoRT_species*i
                    SoR_rxn += SoR_species*i
                GoRT_rxn += GoRT_species*i
        if verbose and not list_calc:
            print(("-"*(10+len(reaction))))
            print(("Total\t\t%.2f\t%.2f\t%.2f" % (HoRT_rxn, SoR_rxn, GoRT_rxn)))
            print(("-"*(10+len(reaction))))
        return GoRT_rxn

    @classmethod
    def from_csv(file_name, verbose = True):
        shomates = []
        with open(file_name, 'r') as shomate_file:
            for line in shomate_file:
                #If the line is not a comment
                if line[0] != '!':
                    data = line.split(',')
                    if verbose:
                        print(("Importing %s" % data[0]))
                    shomates.append(shomate(symbol = data[0], T_low = float(data[1]), T_high = float(data[2]), a = np.array([float(i) for i in data[3:] if i != '' and i != '\n'])))
        return cls(shomates = shomates)

    @classmethod
    def read_shomate(shomate_path, verbose = True):
        shomates = []
        with open(shomate_path, 'r') as shomate_file:
            for line in shomate_file:
                #If the line is not a comment
                if line[0] != '!':
                    data = line.split(',')
                    if verbose:
                        print(("Importing %s" % data[0]))
                    shomates.append(shomate(symbol = data[0], T_low = float(data[1]), T_high = float(data[2]), a = np.array([float(i) for i in data[3:] if i != '' and i != '\n'])))
        return cls(shomates = shomates)
