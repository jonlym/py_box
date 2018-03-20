# -*- coding: utf-8 -*-
"""
Constants module
This module contains universal constants for catalysis research
Created on Wed Nov 23 21:10:42 2016

@author: Jon Lym
"""
import numpy as np
from warnings import warn

def R(units):
    """
    Returns the universal molar gas constant, R.
    Supported unit options: 
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
        'eV/K': 8.6173303e-5
    }
    try:
        return R_dict[units]
    except KeyError:
        raise Exception('Invalid unit: {}'.format(units))
    
def h(units, bar = False):
    """
    Returns the Planck's constant, h. Use bar = True to return h/2π
    Supported unit options: 
    > J s
    > kJ s
    > eV s
    """
    h_dict = {
        'J s': 6.626070040e-34,
        'kJ s': 6.626070040e-37,
        'eV s': 4.135667662e-15
    }

    try:
        h_dict[units]
    except KeyError:
        raise Exception('Invalid unit: {}'.format(units))

    if bar:
        return h_dict[units]/(2.*np.pi)
    else:
        return h_dict[units]

def kb(units):
    """
    Returns the Boltzmann constant.
    Supported unit options: 
    > J/K
    > eV/K
    > cal/K
    
    """
    kb_dict = {
        'J/K': 1.38064852e-23,
        'kJ/K': 1.38064852e-26,        
        'eV/K': 8.6173303e-5,
        'cal/K': 3.2976230e-24,
        'kcal/K': 3.2976230e-27
    }
    try:
        return kb_dict[units]
    except KeyError:
        raise KeyError('Invalid unit: {}'.format(units))

def c(units):
    """
    Returns the speed of light.
    Supported unit options: 
    > m/s
    > cm/s    
    """
    c_dict = {
        'm/s': 299792458.,
        'cm/s': 299792458.e2,
    }
    try:
        return c_dict[units]
    except KeyError:
        raise KeyError('Invalid unit: {}'.format(units))

def m_e(units):
    """
    Returns the mass of an electron.
    Supported unit options: 
    > kg
    > g
    > amu   
    """
    m_e_dict = {
        'kg': 9.10938356e-31,
        'g': 9.10938356e-28,
        'amu': 5.48579909070e-4
    }
    try:
        return m_e_dict[units]
    except KeyError:
        raise KeyError('Invalid unit: {}'.format(units))

def m_p(units):
    """
    Returns the mass of a proton.
    Supported unit options: 
    > kg
    > amu   
    """
    m_p_dict = {
        'kg': 1.6726219e-27,
        'g': 1.6726219e-24,
        'amu': 1.007276466879
    }
    try:
        return m_p_dict[units]
    except KeyError:
        raise KeyError('Invalid unit: {}'.format(units))

def P0(units):
    """
    Returns the reference pressure.
    Supported unit options: 
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
        'bar': 1.01325,
        'atm': 1.,
        'Pa': 101325.,
        'kPa': 101.325,
        'MPa': 0.101325,
        'psi': 14.6959,
        'mmHg': 760.,
        'Torr': 760.,
        None: 1
    }
    try:
        return P0_dict[units]
    except KeyError:
        raise KeyError('Invalid unit: {}'.format(units))

def T0(units):
    """
    Returns room temperature.
    Supported units:
    > K
    > C
    > R
    > F
    """
    T0_dict = {
    'K': 298.15,
    'C': 25.,
    'R': 533.07,
    'F': 73.4
    }
    try:
        return T0_dict[units]
    except KeyError:
        raise KeyError('Invalid unit: {}'.format(units))


Na = 6.02214086e23
e = 1.6021766208e-19

def convert_unit(num = 1., from_ = None, to = None):
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
            |psi        | Pounds per square inch
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
        'hr': 'time',
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
        'torr': 'pressure',
        'psi': 'pressure'
    }
    
    unit_dict = {
        'J': 1.,
        'kJ': 1.e-3,
        'eV': 6.242e+18,
        'cal': 0.239006,
        'kcal': 0.000239006,
        'L atm': 101.33,
        'J/mol': 1.,
        'kJ/mol': 1.e-3,
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
        'cm': 100.,
        'nm': 1.e9,
        'A': 1.e10,
        'm2': 1.,
        'cm2': 1.e4,
        'A2': 1.e20,
        'm3': 1.,
        'cm3': 1.e6,
        'mL': 1.e6,
        'L': 1.e3,
        'kg': 1.,
        'g': 1.e3,
        'amu': 6.022e+26,
        'Pa': 1.,
        'kPa': 1.e-3,
        'MPa': 1.e-6,
        'atm': 9.86923e-6,
        'bar': 1.e-5,
        'mmHg': 0.00750062,
        'torr': 0.00750062,
        'psi': 0.000145038
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