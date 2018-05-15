import collections
import numpy as np
from py_box3.thermo.shomate import Shomate

class Shomates(collections.UserDict):
    """
    A User Dictionary object that stores shomate objects. The key is the symbol name.
    """
    def print_symbols(self):
        """
        Prints a summary of the thermdat list.
        """
        print("Index\tSymbol")
        for i, symbol in enumerate(self.keys()):
            print("{}\t{}".format(i, symbol))

    def property_list(self, property_name):
        """
        Returns the property as a list

        Attributes
        ----------
            property_name - string
            Property to return
        Returns
        -------
            property_list - list
            List of the property requested 
        """
        return [shomate.__dict__[property_name] for shomate in self.values()]

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
