from py_box.cluster_expansion.configuration import Configuration, default_dict
import numpy as np
from warnings import warn

class Configurations(object):
    def __init__(self, configurations = None, info = None):
        self._configurations = []
        if configurations is not None:
            for configuration in configurations:
                self._configurations.append(configuration)

        self.info = info

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
        return np.array([configuration.E_DFT for configuration in self])

    def set_DFT_energies(self, DFT_Es):
        for DFT_E, configuration in zip(DFT_Es, self):
            configuration.E_DFT = DFT_E

    def get_CE_energies(self):
        return np.array([configuration.E_CE for configuration in self])

    def set_CE_energies(self, CE_Es):
        for CE_E, configuration in zip(CE_Es, self):
            configuration.E_CE = CE_E

    def get_copy(self, indices = []):
        """Creates a new Clusters object that will contain the indices requested.
        If no indices are specified, the whole object will be copied."""
        configurations = []
        #If indices not specified, copy whole object
        if len(indices) == 0:
            indices = range(len(self))

        for i in indices:
            configurations.append(copy(self[i]))
        return Configurations(configurations = configurations)

    def get_n_vacancies(self):
        return np.array([configuration.n_vacancies for configuration in self])

    def __len__(self):
        return len(self._configurations)

    def __setitem__(self, index, configuration):
        self._configurations[index] = configuration

    def __getitem__(self, index):
        return self._configurations[index]

    @classmethod
    def from_excel(cls, file_name = 'configurations.xlsx', skiprows = [], info = ''):
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
                                                n_vacancies=configuration['n_vacancies']))
        return cls(configurations = configurations, info = info)
