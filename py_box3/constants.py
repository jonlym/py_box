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
    Universal molar gas constant, R.
    Parameters
        units - str
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
    Returns
        float
            Universal molar gas constant in appropriate units
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
    Planck's constant, h. 
    Parameters
        units - str
            Supported unit options: 
            > J s
            > kJ s
            > eV s
        bar - bool
            If true, returns h/2π       
    Returns
        float
            Planck's constant in appropriate units
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
    Boltzmann constant
    Parameters
        units - str
            Supported unit options: 
            > J/K
            > eV/K
            > cal/K
    Returns
        float
            Boltzmann constant in appropriate units    
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
    Speed of light.
    Parameters
        units - str
            Supported unit options: 
            > m/s
            > cm/s    
    Returns
        float
            Speed of light in appropriate units
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
    Mass of an electron.
    Parameters
        units - str
            Supported unit options: 
            > kg
            > g
            > amu   
    Returns
        float
            Mass of electron in appropriate units
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
    Mass of a proton
    Parameters
        units - str
            Supported unit options: 
            > kg
            > amu   
    Returns
        float
            Mass of proton in appropriate units
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
    Reference pressure.
    Parameters
        units - str
            Supported unit options: 
            > bar
            > atm
            > Pa
            > kPa
            > MPa
            > psi
            > mmHg
            > Torr
    Returns
        float
            Reference pressure in appropriate units
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
    Room temperature.
    Parameters
        units - str
            Supported units:
            > K
            > C
            > R
            > F
    Returns
        float
            Room temperature in appropriate units
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
    Converts units between two unit sets.
    Parameters
        num - float
            Number to convert. This must be specified if converting
            temperatures
        from_ - str
            Units that num is currently in
        to - str
            Units you would like num to be in

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
    Returns
        float
            num in the appropriate units
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


def get_molecular_weight(elements):
    """
    Molecular mass (in g/mol) given the elemental composition.
    Data taken from: https://en.wikipedia.org/wiki/Standard_atomic_weight

    Parameters
        elements - dict
            Elemental composition of species where the keys are the element symbol, atomic number, 
            or element name and the value is the stoichiometric coefficient.

    Returns
        molecular_weight - float
            Molecular weight as float
    """
    molecular_weight = 0.
    for element, coefficient in elements.items():
        molecular_weight += atomic_weight[element] * coefficient
    return molecular_weight

def parse_unit(unit):
    """
    Parses a unit and rewrites it in a standard notation.

    Parameters
        unit - str
            Unit should contain a space between each dimension. e.g. J s
            Indices are allowed but must directly follow the unit. e.g. m2, s-1
            Prefixes are allowed but must directly follow the unit. e.g. kJ
    Returns
        dict
            Standard unit form
    """
    unit_split = unit.split(' ')
    unit_dict = {}
    for unit in unit_split:


        pass

prefixes = {
    'Y': 1.e24,
    'Z': 1.e21,
    'E': 1.e18,
    'P': 1.e15,
    'T': 1.e12,
    'G': 1.e9,
    'M': 1.e6,
    'k': 1.e3,
    'm': 1.e-3,
    'mu': 1.e-9,
    'p': 1.e-12,
    'f': 1.e-15,
    'a': 1.e-18,
    'z': 1.e-21,
    'y': 1.e-24
}

common_units = (
    'C', 'F', 'K', #Temperature
    'l', 'L', 'gal', 'fl', 'pt', #Volume
    'm', 'in', 'ft', 'yd', #Distance
    'g', 'lb', 'amu', #Mass
    's', 'min', 'h', 'day', 'week', 'yr', #Time
    'Hz', #Frequency
    'atm', 'bar', 'Pa', 'psi', 'psia', 'psig', 'torr', 'mmHg', #Pressure
    'J', 'cal', 'btu', #Energy,
    'N', #Force
    'W', #Power
    )


#Dictionary with atomic weights of every element. Elements can be searched
#by atomic number, element symbol or written name in lower case.
atomic_weight = {
    1: 1.008,
    2: 4.002602,
    3: 6.938,
    4: 9.0121831,
    5: 10.806,
    6: 12.0116,
    7: 14.007,
    8: 15.999,
    9: 18.99840316,
    10: 20.1797,
    11: 22.98976928,
    12: 24.305,
    13: 26.9815385,
    14: 28.085,
    15: 30.973762,
    16: 32.06,
    17: 35.45,
    18: 39.948,
    19: 39.0983,
    20: 40.078,
    21: 44.955908,
    22: 47.867,
    23: 50.9415,
    24: 51.9961,
    25: 54.938044,
    26: 55.845,
    27: 58.933194,
    28: 58.6934,
    29: 63.546,
    30: 65.38,
    31: 69.723,
    32: 72.63,
    33: 74.921595,
    34: 78.971,
    35: 79.901,
    36: 83.798,
    37: 85.4678,
    38: 87.62,
    39: 88.90584,
    40: 91.224,
    41: 92.90637,
    42: 95.95,
    43: 98,
    44: 101.07,
    45: 102.9055,
    46: 106.42,
    47: 107.8682,
    48: 112.414,
    49: 114.818,
    50: 118.71,
    51: 121.76,
    52: 127.6,
    53: 126.90447,
    54: 131.293,
    55: 132.905452,
    56: 137.327,
    57: 138.90547,
    58: 140.116,
    59: 140.90766,
    60: 144.242,
    61: 145,
    62: 150.36,
    63: 151.964,
    64: 157.25,
    65: 158.92535,
    66: 162.5,
    67: 164.93033,
    68: 167.259,
    69: 168.93422,
    70: 173.054,
    71: 174.9668,
    72: 178.49,
    73: 180.94788,
    74: 183.84,
    75: 186.207,
    76: 190.23,
    77: 192.217,
    78: 195.084,
    79: 196.966569,
    80: 200.592,
    81: 204.382,
    82: 207.2,
    83: 208.9804,
    84: 209,
    85: 210,
    86: 222,
    87: 223,
    88: 226,
    89: 227,
    90: 232.0377,
    91: 231.03588,
    92: 238.02891,
    93: 237,
    94: 244,
    95: 243,
    96: 247,
    97: 247,
    98: 251,
    99: 252,
    100: 257,
    101: 258,
    102: 259,
    103: 262,
    104: 267,
    105: 268,
    106: 271,
    107: 272,
    108: 270,
    109: 276,
    110: 281,
    111: 280,
    112: 285,
    113: 284,
    114: 289,
    115: 288,
    116: 293,
    118: 294,
    'H': 1.008,
    'He': 4.002602,
    'Li': 6.938,
    'Be': 9.0121831,
    'B': 10.806,
    'C': 12.0116,
    'N': 14.007,
    'O': 15.999,
    'F': 18.99840316,
    'Ne': 20.1797,
    'Na': 22.98976928,
    'Mg': 24.305,
    'Al': 26.9815385,
    'Si': 28.085,
    'P': 30.973762,
    'S': 32.06,
    'Cl': 35.45,
    'Ar': 39.948,
    'K': 39.0983,
    'Ca': 40.078,
    'Sc': 44.955908,
    'Ti': 47.867,
    'V': 50.9415,
    'Cr': 51.9961,
    'Mn': 54.938044,
    'Fe': 55.845,
    'Co': 58.933194,
    'Ni': 58.6934,
    'Cu': 63.546,
    'Zn': 65.38,
    'Ga': 69.723,
    'Ge': 72.63,
    'As': 74.921595,
    'Se': 78.971,
    'Br': 79.901,
    'Kr': 83.798,
    'Rb': 85.4678,
    'Sr': 87.62,
    'Y': 88.90584,
    'Zr': 91.224,
    'Nb': 92.90637,
    'Mo': 95.95,
    'Tc': 98,
    'Ru': 101.07,
    'Rh': 102.9055,
    'Pd': 106.42,
    'Ag': 107.8682,
    'Cd': 112.414,
    'In': 114.818,
    'Sn': 118.71,
    'Sb': 121.76,
    'Te': 127.6,
    'I': 126.90447,
    'Xe': 131.293,
    'Cs': 132.905452,
    'Ba': 137.327,
    'La': 138.90547,
    'Ce': 140.116,
    'Pr': 140.90766,
    'Nd': 144.242,
    'Pm': 145,
    'Sm': 150.36,
    'Eu': 151.964,
    'Gd': 157.25,
    'Tb': 158.92535,
    'Dy': 162.5,
    'Ho': 164.93033,
    'Er': 167.259,
    'Tm': 168.93422,
    'Yb': 173.054,
    'Lu': 174.9668,
    'Hf': 178.49,
    'Ta': 180.94788,
    'W': 183.84,
    'Re': 186.207,
    'Os': 190.23,
    'Ir': 192.217,
    'Pt': 195.084,
    'Au': 196.966569,
    'Hg': 200.592,
    'Tl': 204.382,
    'Pb': 207.2,
    'Bi': 208.9804,
    'Po': 209,
    'At': 210,
    'Rn': 222,
    'Fr': 223,
    'Ra': 226,
    'Ac': 227,
    'Th': 232.0377,
    'Pa': 231.03588,
    'U': 238.02891,
    'Np': 237,
    'Pu': 244,
    'Am': 243,
    'Cm': 247,
    'Bk': 247,
    'Cf': 251,
    'Es': 252,
    'Fm': 257,
    'Md': 258,
    'No': 259,
    'Lr': 262,
    'Rf': 267,
    'Db': 268,
    'Sg': 271,
    'Bh': 272,
    'Hs': 270,
    'Mt': 276,
    'Ds': 281,
    'Rg': 280,
    'Cn': 285,
    'Uut': 284,
    'Fl': 289,
    'Uup': 288,
    'Lv': 293,
    'Uuo': 294,
    'hydrogen': 1.008,
    'helium': 4.002602,
    'lithium': 6.938,
    'beryllium': 9.0121831,
    'boron': 10.806,
    'carbon': 12.0116,
    'nitrogen': 14.007,
    'oxygen': 15.999,
    'fluorine': 18.99840316,
    'neon': 20.1797,
    'sodium': 22.98976928,
    'magnesium': 24.305,
    'aluminium': 26.9815385,
    'silicon': 28.085,
    'phosphorus': 30.973762,
    'sulfur': 32.06,
    'chlorine': 35.45,
    'argon': 39.948,
    'potassium': 39.0983,
    'calcium': 40.078,
    'scandium': 44.955908,
    'titanium': 47.867,
    'vanadium': 50.9415,
    'chromium': 51.9961,
    'manganese': 54.938044,
    'iron': 55.845,
    'cobalt': 58.933194,
    'nickel': 58.6934,
    'copper': 63.546,
    'zinc': 65.38,
    'gallium': 69.723,
    'germanium': 72.63,
    'arsenic': 74.921595,
    'selenium': 78.971,
    'bromine': 79.901,
    'krypton': 83.798,
    'rubidium': 85.4678,
    'strontium': 87.62,
    'yttrium': 88.90584,
    'zirconium': 91.224,
    'niobium': 92.90637,
    'molybdenum': 95.95,
    'technetium': 98,
    'ruthenium': 101.07,
    'rhodium': 102.9055,
    'palladium': 106.42,
    'silver': 107.8682,
    'cadmium': 112.414,
    'indium': 114.818,
    'tin': 118.71,
    'antimony': 121.76,
    'tellurium': 127.6,
    'iodine': 126.90447,
    'xenon': 131.293,
    'cesium': 132.905452,
    'barium': 137.327,
    'lanthanum': 138.90547,
    'cerium': 140.116,
    'praseodymium': 140.90766,
    'neodymium': 144.242,
    'promethium': 145,
    'samarium': 150.36,
    'europium': 151.964,
    'gadolinium': 157.25,
    'terbium': 158.92535,
    'dysprosium': 162.5,
    'holmium': 164.93033,
    'erbium': 167.259,
    'thulium': 168.93422,
    'ytterbium': 173.054,
    'lutetium': 174.9668,
    'hafnium': 178.49,
    'tantalum': 180.94788,
    'tungsten': 183.84,
    'rhenium': 186.207,
    'osmium': 190.23,
    'iridium': 192.217,
    'platinum': 195.084,
    'gold': 196.966569,
    'mercury': 200.592,
    'thallium': 204.382,
    'lead': 207.2,
    'bismuth': 208.9804,
    'polonium': 209,
    'astatine': 210,
    'radon': 222,
    'francium': 223,
    'radium': 226,
    'actinium': 227,
    'thorium': 232.0377,
    'protactinium': 231.03588,
    'uranium': 238.02891,
    'neptunium': 237,
    'plutonium': 244,
    'americium': 243,
    'curium': 247,
    'berkelium': 247,
    'californium': 251,
    'einsteinium': 252,
    'fermium': 257,
    'mendelevium': 258,
    'nobelium': 259,
    'lawrencium': 262,
    'rutherfordium': 267,
    'dubnium': 268,
    'seaborgium': 271,
    'bohrium': 272,
    'hassium': 270,
    'meitnerium': 276,
    'darmstadtium': 281,
    'roentgenium': 280,
    'copernicium': 285,
    'ununtrium': 284,
    'flerovium': 289,
    'ununpentium': 288,
    'livermorium': 293,
    'ununoctium': 294,
}