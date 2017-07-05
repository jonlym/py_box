import numpy as np
from ase import Atom as Atom_ase

class Atom(Atom_ase):
    """"
    Modified Atom class from ASE.
    https://wiki.fysik.dtu.dk/ase/ase/atom.html

    Added attributes:
    bader_charge: Float
        Bader charge of atom in e.
    dos_occupancy: Dict
        Electron occupancy with the orbital as the key and the value as the value in e.
    cn: Int
        Coordination number.
    gcn: Float
        Generalized coordination number
    neighbors: Set
        Indices of coordinated atoms
    """

    def __init__(self, symbol = 'X', position = (0, 0, 0), tag = None, momentum = None, mass = None, magmom = None,
                 charge = None, atoms = None, index = None, bader_charge = None, dos_occupancy = None, cn = None,
                 gcn = None, neighbors = None):
        Atom_ase.__init__(self, symbol=symbol, position = position, tag = tag, momentum = momentum, mass = mass, magmom = magmom, charge = charge, atoms = atoms, index = index)
        self.bader_charge = bader_charge
        self.dos_occupancy = dos_occupancy
        self.cn = cn
        self.gcn = gcn
        self.neighbors = neighbors

    def dos_occupancy_total(self):
        return np.sum(self.dos_occupancy.values())
