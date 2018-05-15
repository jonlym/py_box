# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 21:10:42 2016

@author: Jon Lym
"""

import itertools
import numpy as np
from datetime import datetime
import pickle

hex_to_bin_dict = {'0': '0000',
                   '1': '0001',
                   '2': '0010',
                   '3': '0011',
                   '4': '0100',
                   '5': '0101',
                   '6': '0110',
                   '7': '0111',
                   '8': '1000',
                   '9': '1001',
                   'a': '1010',
                   'b': '1011',
                   'c': '1100',
                   'd': '1101',
                   'e': '1110',
                   'f': '1111'}

def convert_hex_to_bin(hex_string, n = None):
    #Add initial zeros
    if n is None:
        out = ''
    else:
        out = '0'*(n - 4*len(hex_string))

    hex_string = hex_string.replace('0x', '').replace('L', '')
    for val in hex_string:
        out = '{}{}'.format(out, hex_to_bin_dict[val])
    return out

def basen_to_base10(num, n):
    out = 0
    length = len(num)
    for i, val in enumerate(num):
        out += int(val) * (n**(length - i - 1))
    return out

def base10_to_basen(num, n, width):
    out = np.zeros(shape = (width, ))
    if num == 0:
        return out
    i = 0
    while num:
        num, r = divmod(num, n)
        out[i] = r
        i += 1
    return out[::-1]

def any_alpha(string):
    """
    Returns True if any alphabetic characters are in the string. False otherwise.
    """
    for character in string:
        if character.isalpha():
            return True
    else:
        return False

def get_unique_list(data):
    """
    Given a list, returns a unique list.
    """
    keys = {}
    for item in data:
        keys[item] = 1
    return list(keys.keys())

def plot_parity(x, y, decimals = 2):
    """
    Plots a party plot given two vectors of data.
    """
    from matplotlib import pyplot as plt
    data = np.concatenate((x, y), axis = 0)
    min_data = np.round(min(data*10.**decimals))/10.**decimals
    max_data = np.round(max(data*10.**decimals))/10.**decimals
    fig = plt.figure()
    plt.plot(np.array([min_data, max_data]), np.array([min_data, max_data]), 'k-', x, y, 'bo')
    axes = plt.gca()
    axes.set_ylim([min_data, max_data])
    axes.set_xlim([min_data, max_data])
    return (fig, axes)

def get_time():
    """
    Returns the time.
    Example: 2017-12-09 23:02:16.459000
    """
    return str(datetime.now())

def get_null(mat, rtol = 1.e-5):
    """
    Returns the nullspace of a 2D matrix, mat
    """
    u, s, v = np.linalg.svd(mat)
    rank = (s > rtol*s[0]).sum()
    return v[rank:].T.copy()

def get_MSE(xs_data, xs_fit):
    """
    Returns the root mean squared error given two vectors of equal lengths
    """
    return np.mean([(x_data-x_fit)**2 for x_data, x_fit in zip(xs_data, xs_fit)])

def get_RMSE(xs_data, xs_fit):
    """
    Returns the root mean squared error given two vectors of equal lengths
    """
    return np.sqrt(get_MSE(xs_data = xs_data, xs_fit = xs_fit))

def spherical_to_xyz(r = 1., theta = 0., phi = 0., degrees = True):
    """
    Converts spherical coordinates to Cartesian coordinates. Angles are in degrees by default. Set degrees to False
    to use radians
    """
    if degrees:
        theta = np.radians(theta)
        phi = np.radians(phi)
    x = r * np.cos(theta) * np.sin(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(phi)
    return np.array([x, y, z])

def get_n_blanks(n):
    return ' '*n

def dict_products(dicts):
    return (dict(zip(dicts, x)) for x in itertools.product(*dicts.values()))

def interpolate(x_low, x_high, y_low, y_high, x):
    if x == x_low:
        return y_low
    elif x == x_high:
        return y_high

    m = (y_high - y_low)/(x_high - x_low)
    b = y_high - m * x_high
    return m * x + b

def saveplt(ax, filename):
    with open(filename, 'wb') as f_ptr:
        pickle.dump(ax, f_ptr)

def loadplt(filename, show = False):
    from matplotlib import pyplot as plt

    with open(filename, 'rb') as f_ptr:
        ax = pickle.load(f_ptr)
    if show:
        plt.show()
    else:
        return ax

def get_molecular_weight(elements):
    molecular_weight = 0.
    for element, coefficient in elements.items():
        molecular_weight += atomic_weight[element] * coefficient
    return molecular_weight

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