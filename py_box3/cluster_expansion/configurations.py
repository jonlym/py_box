from py_box3.cluster_expansion.configuration import Configuration, default_dict
from py_box3.cluster_expansion.In2O3 import get_sigma_from_sites
import numpy as np
import os
from copy import copy
from ase.io import read
from warnings import warn
import platform
from itertools import permutations
from py_box3.ase import DFT_E_gas
try:
    import networkx as nx
except:
    pass
    
class Configurations(object):
    def __init__(self, configurations = None, info = None, del_E_gas = 0.):
        self._configurations = list(configurations)
        self.info = info
        self.del_E_gas = del_E_gas

    def append(self, configuration):
        self._configurations.append(configuration)

    def extend(self, configurations):
        self._configurations.extend(configurations)

    def index(self, name):
        for i, configuration in enumerate(self._configurations):
            if configuration.name == name:
                return i
        else:
            return None

    def remove(self, name = None, index = None):
        if name is not None and index is not None:
            raise Exception('Both name and index cannot be specified.')
        elif name is None and index is None:
            raise Exception('Either name or index must be specified.')
        elif name is not None:
            index = self.index(name)
        del self._configurations[index]

    def get_DFT_energies(self):
        """
        Returns the DFT energies as a numpy array.
        """
        return np.array([configuration.E_DFT for configuration in self])

    def set_DFT_energies(self, DFT_Es):
        """
        Assigns DFT energies based on iterable object DFT_Es. Assigns by element.
        """
        for DFT_E, configuration in zip(DFT_Es, self):
            configuration.E_DFT = DFT_E

    def get_CE_energies(self):
        """
        Returns the cluster expansion energies as a numpy array.
        """
        return np.array([configuration.E_CE for configuration in self])

    def set_CE_energies(self, CE_Es):
        """
        Assigns Cluster Expansion energies based on the iterable object CE_Es. Assigns by element.
        """
        for CE_E, configuration in zip(CE_Es, self):
            configuration.E_CE = CE_E

    def get_copy(self, indices = []):
        """Creates a new Clusters object that will contain the indices requested.
        If no indices are specified, the whole object will be copied."""
        configurations = []
        #If indices not specified, copy whole object
        if len(indices) == 0:
            indices = list(range(len(self)))

        for i in indices:
            configurations.append(copy(self[i]))
        return Configurations(configurations = configurations)

    def write_to_excel(self, file_name):
        """
        Writes the configurations object to a spreadsheet.
        """
        import pandas as pd
        writer = pd.ExcelWriter(file_name)
        configs_data_frame = self._get_data_frame()
        column_order = list(configs_data_frame.columns.values)
        #Inserting name first, sigma second
        column_order.remove('name')
        column_order.insert(0, 'name')
        column_order.remove('sigma')
        column_order.insert(1, 'sigma')

        #Writing to spreadsheet
        configs_data_frame.to_excel(writer, 'Sheet1', columns = column_order, index = False)
        writer.save()

    def _get_data_frame(self):
        """
        Converts the configuration to a Pandas dataframe
        """
        import pandas as pd
        data_dict = {}
        for key in list(self[0].__dict__.keys()):
            print(key)
            if key == 'sigma':
                #Convert sigma to string separated by ', '
                data_dict[key] = [', '.join(str(x) for x in config.__dict__[key]) for config in self]
            else:
                data_dict[key] = [config.__dict__[key] for config in self]
        return pd.DataFrame(data=data_dict)

    def get_n_nodes(self):
        """
        Returns a numpy array of the number of vacancies for each configurations.
        """
        return np.array([configuration.n_nodes for configuration in self])

    def find_n_max(self):
        """
        Finds the configuration with the largest number of vacancies.
        """
        n_max = max(self.get_n_nodes())
        self.n_max = n_max

    def get_n_max(self):
        return self.n_max

    def find_E_max(self, update = True):
        if update:
            self.find_n_max()

        for configuration in self:
            if configuration.n_nodes == self.n_max:
                self.E_max = configuration.E_DFT

    def find_E_min(self, update = True):
        for configuration in self:
            if configuration.n_nodes == 0:
                self.E_min = configuration.E_DFT

    def _update(self):
        self.find_n_max()
        self.find_E_max()
        self.find_E_min()

    def calc_E_fit(self, update = True):
        """
        Calculates the energy used in cluster expansions based on the formula:
        E = E_DFT - (1- n/n_max)E_DFT_0 - (n/n_max)E_DFT_max

        Assumes that the maximum coverage is already present in the Configurations object.
        """
        if update:
            self._update()

        for i, configuration in enumerate(self):
            n = configuration.n_nodes
            self[i].E_fit = convert_DFT_to_CE(E = configuration.E_DFT, E_min = self.E_min, E_max = self.E_max, n = n, n_max = self.n_max)

    def get_E_fit(self, update = True):
        """
        Returns the energy used for cluster expansions in a numpy array.
        Update
         """
        if update:
            self.calc_E_fit(update = update)
        return np.array([configuration.E_fit for configuration in self])

    def get_bit_strings(self):
        X = np.array([configuration.sigma for configuration in self])
        X[X == -1] = 0
        return X

    def make_network(self, calculate_energy = False, use_CE = True, update = True):
        """
        Creates a graph based map that connects configurations with subsequent vacancies
        """
        if update:
            self._update()
        self.network = nx.DiGraph()
        self.network.add_nodes_from([configuration.name for configuration in self])
        for pair in permutations(self, 2):
            config1 = pair[0]
            config2 = pair[1]

            if config1.n_nodes == (config2.n_nodes - 1):
                n_change = np.sum([1 for s1, s2 in zip(config1.sigma, config2.sigma) if s1 != s2])
                #If only one change has been made
                if n_change == 1:
                    if calculate_energy:
                        #Add vacancy formation energy to graph
                        if use_CE:
                            #Use cluster Expansion energy
                            E1 = convert_CE_to_DFT(n = config1.n_nodes, n_max = self.n_max, E = config1.E_CE, E_min = self.E_min, E_max = self.E_max)
                            E2 = convert_CE_to_DFT(n = config2.n_nodes, n_max = self.n_max, E = config2.E_CE, E_min = self.E_min, E_max = self.E_max)
                        else:
                            #Use DFT energy
                            E1 = config1.E_DFT
                            E2 = config2.E_DFT
                        self.network.add_edge(config1.name, config2.name, del_E = E2 - E1 + self.del_E_gas)
                    else:
                        self.network.add_edge(config1.name, config2.name)                    
    def __len__(self):
        return len(self._configurations)

    def __setitem__(self, index, configuration):
        self._configurations[index] = configuration

    def __getitem__(self, index):
        return self._configurations[index]

    @classmethod
    def from_excel(cls, file_name = 'configurations.xlsx', skiprows = [], info = '', del_E_gas = 0.):
        import pandas as pd

        if info == '':
            info = 'Data generated using class method, from_excel(), using file: {}. Rows skipped: {}'.format(file_name, skiprows)
        configuration_data = pd.read_excel(file_name, skiprows = skiprows, header=0)
        configurations = []
        for i, configuration in configuration_data.iterrows():
            for col in configuration_data.columns:
                #Process blank records
                if type(configuration[col]) is float and np.isnan(configuration[col]):
                    configuration.set_value(col, None)
                sigma = [int(x) for x in configuration['sigma'].split(',')]
            configurations.append(Configuration(name = configuration['name'],
                                                sigma = sigma,
                                                E_CE = configuration['E_CE'],
                                                E_DFT = configuration['E_DFT'],
                                                n_nodes=configuration['n_nodes']))
        return cls(configurations = configurations, info = info, del_E_gas = del_E_gas)

    @classmethod
    def from_vasp(cls, path = '.', info = '', del_E_gas = 0.):
        configurations = []
        for root, folders, files in os.walk(path):
            for file in files:
                if 'OUTCAR' in file:
                    try:
                        atoms = read(os.path.join(root, file))
                    except:
                        continue
                    else:
                        full_path = os.path.join(root, file)
                        if 'windows' in platform.system().lower():
                            name = full_path.split('\\')[-2]
                        else:
                            name = full_path.split('/')[-2]
                        E_DFT = atoms.get_potential_energy()
                        if 'clean' in name:
                            n_nodes = 0
                        else:
                            n_nodes = len(name.split('_'))
                        sigma = get_sigma_from_sites(name)

                        configurations.append(Configuration(name = name,
                                                            sigma = sigma,
                                                            E_DFT = E_DFT,
                                                            n_nodes = n_nodes,
                                                            atoms = atoms))
        return cls(configurations = configurations, info = info, del_E_gas = del_E_gas)

def convert_CE_to_DFT(n, n_max, E, E_min, E_max):
    return E + (1. - float(n)/float(n_max)) * E_min + float(n)/float(n_max) * E_max

def convert_DFT_to_CE(n, n_max, E, E_min, E_max):
    #return E - (1. - float(n)/float(n_max)) * E_min - float(n)/float(n_max) * E_max    
    return (E - E_min + n*(DFT_E_gas['H2O'] - DFT_E_gas['H2']))/n_max    