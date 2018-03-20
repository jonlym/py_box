from ase.thermochemistry import IdealGasThermo
import numpy as np
from py_box3 import constants as c

def get_chemical_potential(atoms,
                           freq,
                           geometry,
                           potentialenergy,
                           symmetrynumber,
                           spin,
                           T,
                           nasa_liq = None,
                           H_gas = None,
                           S_gas = None,
                           x = None,
                           gas_phase = False,
                           verbose = False):
    """
    Calculates the chemical potential of the toluene-phase species
    atoms - ASE Atoms object
    freq - Numpy array of floats
        Frequencies of the species in 1/cm
    geometry - string
        Geometry of the species (monatomic, linear, nonlinear)
    potentialenergy - float
        DFT Energy of the species in eV
    symmetrynumber - float
        Symmetry number of the species
    spin - int
        Number of unpaired electrons in the species
    nasa_liq - NASA object
        NASA polynomial of the species in the liquid phase
    H_gas - float
        Enthalpy of the gas in kJ/mol/K
    S_gas - float
        Entropy of the gas in J/mol/K
    T - float
        Temperature to perform analysis
    x - float
        Mole fraction of species in liquid
    """

    #DFT Delta G
    DFT_thermo = IdealGasThermo(vib_energies=freq*c.h('eV s')*c.c('cm/s'),
                                geometry=geometry,
                                potentialenergy=potentialenergy,
                                natoms = len(atoms),
                                atoms = atoms,
                                symmetrynumber=symmetrynumber,
                                spin = spin)
    G_DFT = DFT_thermo.get_gibbs_energy(temperature = 298., pressure = 1.e5, verbose = verbose)

    if gas_phase:
        return np.array([DFT_thermo.get_gibbs_energy(temperature = T_i, pressure = 1.e5, verbose = verbose) for T_i in T])
    else:
        #ASPEN Liquid Data
        G_liq = nasa_liq.get_GoRT(T = T, verbose = False) * c.kb('eV/K') * T

        #NIST Gas Phase Data
        G_gas = H_gas * c.convert_unit(from_ = 'kJ/mol', to = 'eV/molecule') - T * S_gas * c.convert_unit(from_ = 'J/mol', to = 'eV/molecule')

        return G_DFT + G_liq - G_gas + c.kb('eV/K') * T * np.log(x)
