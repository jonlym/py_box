from py_box.cluster_expansion.configuration import Configuration, default_dict
from py_box.cluster_expansion.In2O3 import get_sigma_from_sites
import numpy as np
import os
from copy import copy
from ase.io import read
from warnings import warn
import platform

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

    def write_to_excel(self, file_name):
        """
        Writes the configurations object to a spreadsheet.
        :param file_name: Name of the spreadsheet.
        :return:
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
        :return:
        """
        import pandas as pd
        data_dict = {}
        for key in self[0].__dict__.iterkeys():
            print key
            if key == 'sigma':
                #Convert sigma to string separated by ', '
                data_dict[key] = [', '.join(str(x) for x in config.__dict__[key]) for config in self]
            else:
                data_dict[key] = [config.__dict__[key] for config in self]
        return pd.DataFrame(data=data_dict)

    def get_n_vacancies(self):
        return np.array([configuration.n_vacancies for configuration in self])

    def set_E_fit(self):
        n_max = 0
        for i, configuration in enumerate(self):
            print configuration.name
            if configuration.n_vacancies == 0:
                E0 = configuration.E_DFT

            if configuration.n_vacancies > n_max:
                E_max = configuration.E_DFT
                n_max = configuration.n_vacancies

        for i, configuration in enumerate(self):
            n = configuration.n_vacancies
            self[i].E_fit = configuration.E_DFT - (1.-float(n)/float(n_max))*E0 - float(n)/float(n_max)*E_max

    def get_E_fit(self):
        return np.array([configuration.E_fit for configuration in self])

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

    @classmethod
    def from_vasp(cls, path = '.', info = ''):
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
                            n_vacancies = 0
                        else:
                            n_vacancies = len(name.split('_'))
                        sigma = get_sigma_from_sites(name)

                        configurations.append(Configuration(name = name,
                                                            sigma = sigma,
                                                            E_DFT = E_DFT,
                                                            n_vacancies = n_vacancies))
        return cls(configurations = configurations, info = info)
