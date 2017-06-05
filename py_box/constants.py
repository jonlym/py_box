# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 21:10:42 2016

@author: Jon Lym
"""
import numpy as np
from warnings import warn

def R(units = None):
    """
    Returns the universal molar gas constant. Default unit set is J/mol/K
    Acceptable unit options: 
    > J/mol/K 
    > kJ/mol/K
    > L kPa/mol/K
    > cm3 kPa/mol/K
    > m3 Pa/mol/K
    > cm3 MPa/mol/K
    > m3 bar/mol/K
    > L bar/mol/K
    > L torr/mol/K
    > cal/mol/K
    > kcal/mol/K
    > L atm/mol/K
    > cm3 atm/mol/K
    > eV/K (This option implies per molecule. i.e. Uses the value of kb)
    """
    R_dict = {
        'J/mol/K': 8.3144598,
        'kJ/mol/K': 8.3144598e-3,
        'L kPa/mol/K': 8.3144598,
        'cm3 kPa/mol/K': 8.3144598e3,
        'm3 Pa/mol/K': 8.3144598,
        'cm3 MPa/mol/K': 8.3144598,
        'm3 bar/mol/K': 8.3144598e-5,
        'L bar/mol/K': 8.3144598e-2,
        'L torr/mol/K': 62.363577,
        'cal/mol/K': 1.9872036,
        'kcal/mol/K': 1.9872036e-3,
        'L atm/mol/K': 0.082057338,
        'cm3 atm/mol/K': 82.057338,
        'eV/K': 8.6173303e-5,
        None: 8.3144598
    }
    if units is None:
        warn("Unit type not for R specified. Using default unit set J/mol/K")
    if R_dict.get(units) == None:
        warn("Unit type %r for R not found. Using default unit set J/mol/K" % units)
        units = 'J/mol/K'
    return R_dict[units]
    
def h(units = None, bar = False):
    """
    Returns the Planck's constant. Default unit set is J s. Use bar = True to return h/2π
    Acceptable unit options: 
    > J s
    > kJ s
    > eV s
    """
    h_dict = {
        'J s': 6.626070040e-34,
        'kJ s': 6.626070040e-37,
        'eV s': 4.135667662e-15,
        None: 6.626070040e-34
    }
    if units is None:
        warn("Unit type for h not specified. Using default unit set J s")
    if h_dict.get(units) is None:
        warn("Unit type %r for h not found. Using default unit set J s" % units)
        units = 'J s'
    if bar:
        return h_dict[units]/(2*np.pi)
    else:
        return h_dict[units]

def kb(units = None):
    """
    Returns the Boltzmann constant. Default unit set is J/K
    Acceptable unit options: 
    > J/K
    > eV/K
    > cal/K
    
    """
    kb_dict = {
        'J/K': 1.38064852e-23,
        'kJ/K': 1.38064852e-26,        
        'eV/K': 8.6173303e-5,
        'cal/K': 3.2976230e-24,
        'kcal/K': 3.2976230e-27,
        None: 6.626070040e-34
    }
    if units is None:
        warn("Unit type not kb specified. Using default unit set J/K")
    if kb_dict.get(units) is None:
        warn("Unit type %r for kb not found. Using default unit set J/K" % units)
        units = 'J/K'
    return kb_dict[units]

def c(units = None):
    """
    Returns the speed of light. Default unit set is m/s
    Acceptable unit options: 
    > m/s
    > cm/s    
    """
    c_dict = {
        'm/s': 299792458,
        'cm/s': 299792458e2,
        None: 299792458
    }
    if units is None:
        warn("Unit type not c specified. Using default unit set m/s")
    if c_dict.get(units) is None:
        warn("Unit type %r for c not found. Using default unit set m/s" % units)
        units = 'm/s'
    return c_dict[units]

def m_e(units = None):
    """
    Returns the mass of an electron. Default unit set is kg
    Acceptable unit options: 
    > kg
    > amu   
    """
    m_e_dict = {
        'kg': 9.10938356e-31,
        'amu': 5.48579909070e-4,
        None: 9.10938356e-31
    }
    if units is None:
        warn("Unit type not m_e specified. Using default unit set kg")
    if m_e_dict.get(units) is None:
        warn("Unit type %r for m_e not found. Using default unit set kg" % units)
        units = 'kg'
    return m_e_dict[units]

def m_p(units = None):
    """
    Returns the mass of a proton. Default unit set is kg
    Acceptable unit options: 
    > kg
    > amu   
    """
    m_p_dict = {
        'kg': 1.6726219e-27,
        'amu': 1.007276466879,
        None: 1.6726219e-27
    }
    if units is None:
        warn("Unit type not m_p specified. Using default unit set kg")
    if m_p_dict.get(units) is None:
        warn("Unit type %r for m_p not found. Using default unit set kg" % units)
        units = 'kg'
    return m_p_dict[units]

def P0(units = None):
    """
    Returns the reference pressure. Default unit set is bar
    Acceptable unit options: 
    > bar
    > atm
    > Pa
    > kPa
    > MPa
    > psi
    > mmHg
    > Torr
    """
    P0_dict = {
        'bar': 1,
        'atm': 0.987,
        'Pa': 1e5,
        'kPa': 100.,
        'MPa': 0.1,
        'psi': 14.5038,
        'mmHg': 750.06,
        'Torr': 750.06,
        None: 1
    }
    if units is None:
        warn("Unit type not c specified. Using default unit set bar")
    if P0_dict.get(units) is None:
        warn("Warning: Unit type %r for c not found. Using default unit set m/s" % units)
        units = 'bar'
    return P0_dict[units]


Na = 6.02214086e23
e = 1.6021766208e-19
T0 = 298.

def convert_unit(num = 1, from_ = None, to = None):
    """
    Converts num with units 'from_' to num with units 'to'.
    Supported Units:
    -----------------------------------    
    Type    | Symbol    | Unit
    -----------------------------------
    Energy  |J          | Joules
            |kJ         | KiloJoules
            |eV         | Electron Volts
            |cal        | Calories
            |kcal       | Kilocalories
            |L atm      | Litre atmospheres
    -----------------------------------
    Energy/ |J/mol      |Joules per mol
    Amount  |kJ/mol     |KiloJoules per mol
            |cal/mol    |Calories per mole
            |kcal/mol   |Kilocalories per mole
            |eV/molecule|eV per molecule
    -----------------------------------
    Time    |s          | Seconds 
            |min        | Minutes
            |hr         | Hours
    -----------------------------------
    Amount  |molecule   | Molecule
            |mol        | Moles
    -----------------------------------
    Temp    |C          | Celcius
            |K          | Kelvin
    -----------------------------------
    Length  |m          | Meter
            |cm         | Centimeter
            |nm         | Nanometer           
            |A          | Anstroms
    -----------------------------------
    Area    |m2         | Meters Squared
            |cm2        | Centimeter Squared
            |A2         | Anstroms Squared
    -----------------------------------
    Volume  |m3         | Meters Cubed
            |cm3        | Centimeters Cubed
            |mL         | Milliliters
            |L          | Liters
    -----------------------------------
    Mass    |kg         | Kilograms
            |g          | Grams
            |amu        | Atomic mass units
    -----------------------------------
    Pressure|Pa         | Pascals
            |kPa        | KiloPascals
            |MPa        | MegaPascals
            |atm        | Atmospheres
            |bar        | Bars
            |mmHg       | Millimeters of Mercury
    """
    
    type_dict = {
        'J': 'energy',
        'kJ': 'energy',
        'eV': 'energy',
        'cal': 'energy',
        'kcal': 'energy',
        'L atm': 'energy',
        'J/mol': 'energy/amount',
        'kJ/mol': 'energy/amount',
        'cal/mol': 'energy/amount',
        'kcal/mol': 'energy/amount',
        'eV/molecule': 'energy/amount',
        's': 'time',
        'min': 'time',
        'hr': 'hour',
        'molecule': 'amount',
        'mol': 'amount',
        'C': 'temp',
        'K': 'temp',
        'm': 'length',
        'cm': 'length',
        'nm': 'length',
        'A': 'length',
        'm2': 'area',
        'cm2': 'area',
        'A2': 'area',
        'm3': 'volume',
        'cm3': 'volume',
        'mL': 'volume',        
        'L': 'volume',
        'kg': 'mass',
        'g': 'mass',
        'amu': 'mass',
        'Pa': 'pressure',
        'kPa': 'pressure',
        'MPa': 'pressure',
        'atm': 'pressure',
        'bar': 'pressure',
        'mmHg': 'pressure',
        'torr': 'pressure'
    }
    
    unit_dict = {
        'J': 1.,
        'kJ': 1e-3,
        'eV': 6.242e+18,
        'cal': 0.239006,
        'kcal': 0.000239006,
        'L atm': 101.33,
        'J/mol': 1.,
        'kJ/mol': 1e-3,
        'cal/mol': 0.239006,
        'kcal/mol': 0.000239006,
        'eV/molecule': 6.242e+18/6.02214086e23,
        's': 1.,
        'min': 1./60.,
        'hr': 1./3600.,
        'mol': 1.,
        'molecule': 6.02214086e23,
        'C': 0.,
        'K': 273.15,
        'm': 1.,
        'cm': 100,
        'nm': 1e9,
        'A': 1e10,
        'm2': 1.,
        'cm2': 1e4,
        'A2': 1e20,
        'm3': 1.,
        'cm3': 1e6,
        'mL': 1e6,
        'L': 1e3,
        'kg': 1.,
        'g': 1e3,
        'amu': 6.022e+26,
        'Pa': 1.,
        'kPa': 1e-3,
        'MPa': 1e-6,
        'atm': 9.86923e-6,
        'bar': 1e-5,
        'mmHg': 0.00750062,
        'torr': 0.00750062
    }

    #Check if the entry exists    
    if type_dict.get(from_) is None:
        raise ValueError("%r not a supported unit." % from_)
    if type_dict.get(to) is None:
        raise ValueError("%r not a supported unit." % to)
    #Check that the unit types are the same
    from_type = type_dict[from_]
    to_type = type_dict[to]
    if from_type != to_type:
        raise ValueError("%r [Type %r] not compatible with %r [Type %r]" % (from_, from_type, to, to_type))
    elif from_type == 'temp':
        return num + unit_dict[to] - unit_dict[from_]
    else:
        return num * unit_dict[to] / unit_dict[from_]        
    
    