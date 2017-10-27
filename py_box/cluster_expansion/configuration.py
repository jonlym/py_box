from ase.io import read
import numpy as np
import os

class Configuration(object):
    """
    Holds the configuration.
        Attributes:
            sigma - list of int
                Vector that holds the state of the site. Usually either +1 or -1.
            E_DFT - float
                Energy of the configuration obtained by DFT
            E_CE - float
                Energy of the configuration obtained by cluster expansion
            n_nodes - int
                Number of nodes in cell
    """

    def __init__(self, name = None, sigma = None, E_DFT = None, E_fit = None, E_CE = None, n_nodes = None):
        self.name = name
        self.sigma = sigma
        self.E_DFT = E_DFT
        self.E_fit = E_fit
        self.E_CE = E_CE
        self.n_nodes = n_nodes

    @classmethod
    def from_vasp(cls, name = None, sigma = None, n_nodes = None, path = './'):
        atoms = read(os.path.join(path, 'OUTCAR'))
        try:
            E_DFT = atoms.get_potential_energy()
        except:
            E_DFT = None
        return cls(name = name, sigma = sigma, E_DFT = E_DFT, n_nodes = n_nodes)

default_dict = {'name': str,
                'sigma': list,
                'E_DFT': float,
                'E_CE': float,
                'E_DFT_Raw': float,
                'n_nodes': int}
