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
    """
    Converts hexidemial string to binary string

    Parameters
    ----------
        hex_string - string
            String of the hexidecimal number to convert. 
        n - int
            Length of binary string returned
    Returns
    -------
        bin_string - string
            Binary string
    """    
    #Add initial zeros
    if n is None:
        bin_string = ''
    else:
        bin_string = '0'*(n - 4*len(hex_string))

    hex_string = hex_string.replace('0x', '').replace('L', '')
    for val in hex_string:
        bin_string = '{}{}'.format(bin_string, hex_to_bin_dict[val.lower()])
    return bin_string

def basen_to_base10(num, n):
    """
    Converts a number from base n to base 10.

    Parameters
    ----------
        num - int, string, list, or ndarray
            Base n number to convert
        n - int
            Base of number
    Returns
    -------
        num_10 - int
            Base 10 representation of num_n
    """

    if isinstance(num, int) or isinstance(num, str):
        num = np.array([int(x) for x in str(num)])

    num_10 = 0
    for i, val in enumerate(num):
        num_10 += int(val) * (n**(len(num) - i - 1))
    return num_10

def base10_to_basen(num, n, width = None):
    """
    Converts a number from base 10 to base n.

    Parameters
    ----------
        num - int
            Base 10 number to convert
        n - int
            Base to convert the number
        width - int
            Width of the number returned

    Returns
    -------
        num_n - (N,) ndarray
            Base n number
    """
    num_n_rev = []
    while num >= 0:
        num, r = divmod(num, n)
        num_n_rev.append(r)

    #Pad with zeros
    if width is not None:
        num_n_rev = num_n_rev + [0] * (len(num_n_rev) - width)
    num_n = np.array(num_n_rev)[::-1]
    return num_n

def any_alpha(string):
    """
    Returns True if any alphabetic characters are in the string. False otherwise.

    Parameters
    ----------
        string - string
            String to test for alphabetic characters
    Returns
    -------
        alpha_present - bool
            True if alphabetic characters are present. False otherwise
    """
    for character in string:
        if character.isalpha():
            return True
    else:
        return False

def get_unique_list(data):
    """
    Given a list, returns a list without duplicates.

    Parameters
    ----------
        data - list
            Data with duplicates
    Returns
    -------
        unique_data - list
            Data without duplicates
    """
    keys = {}
    for item in data:
        keys[item] = 1
    return list(keys.keys())

def plot_parity(x, y, decimals = 2, min_val = None, max_val = None):
    """
    Plots a party plot given two vectors of data.

    Parameters
    ----------
        x - (N,) ndarray
            Array of values to plot on the x axis
        y - (N,) ndarray
            Array of values to plot on the y axis
        decimals - int
            Number of decimal points to include on plot
        min_val - float
            Manually override lower limit for parity plot
        max_val - float
            Manually override upper limit for parity plot
    Returns
    -------
        fig - matplotlib.figure.Figure object
        axes - matplotlib.axes.Axes object
    """
    from matplotlib import pyplot as plt

    data = np.concatenate((x, y), axis = 0)
    if min_val is None:
        min_val = np.round(min(data*10.**decimals))/10.**decimals
    if max_val is None:
        max_val = np.round(max(data*10.**decimals))/10.**decimals

    fig = plt.figure()
    plt.plot(np.array([min_val, max_val]), np.array([min_val, max_val]), 'k-', x, y, 'bo')
    axes = plt.gca()
    axes.set_xlim([min_val, max_val])
    axes.set_ylim([min_val, max_val])
    return (fig, axes)

def get_time():
    """
    Returns the time.
    Example: 2017-12-09 23:02:16.459000

    Returns
    -------
        time - str
            Current time
    """
    return str(datetime.now())

def get_null(mat, rtol = 1.e-5):
    """
    Returns the nullspace of a 2D matrix, mat

    Parameters
    ----------
        mat - (M, N) ndarray
            Matrix to find the nullspace
        rtol - float
            Tolerance
    Returns
    -------
        null_space - (N - rank(mat)) ndarray
            Null space of mat
    """
    u, s, v = np.linalg.svd(mat)
    rank = (s > rtol*s[0]).sum()
    return v[rank:].T.copy()

def get_MSE(x, y):
    """
    Returns the root mean squared error given two vectors of equal lengths
    
    Parameters
    ----------
        x - (N,) ndarray
        y - (N,) ndarray

    Returns
        MSE - float
            Mean squared error
    """
    if len(x) != len(y):
        raise ValueError('x and y are not the same length.')
    return np.mean([(x_ - y_)**2 for x_, y_ in zip(x, y)])

def get_RMSE(x, y):
    """
    Returns the root mean squared error given two vectors of equal lengths
    Parameters
    ----------
        x - (N,) ndarray
        y - (N,) ndarray

    Returns
        RMSE - float
            Root mean squared error
    """
    return np.sqrt(get_MSE(x = x, y = y))

def spherical_to_xyz(r = 1., theta = 0., phi = 0., degrees = True):
    """
    Converts spherical coordinates to Cartesian coordinates. Angles are in degrees by default. Set degrees to False
    to use radians

    Parameters
    ----------
        r - float
            Radius
        theta - float
            Azimuthal angle (i.e. Angle between the x and y axes)
        phi - float
            Polar angle (i.e. Angle between the z and the x-y plane)
        degrees - bool
            Whether to use degrees or not. True for degrees, False for radians
    Returns
    -------
        xyz_coordinates - (3,) ndarray
            Cartesian (x, y, z) coordinates
    """
    if degrees:
        theta = np.radians(theta)
        phi = np.radians(phi)
    x = r * np.cos(theta) * np.sin(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(phi)
    return np.array([x, y, z])

def get_n_blanks(n):
    """
    Returns a string made of n blanks

    Parameters
    ----------
        n - int
            Number of blank spaces
    Returns
    -------
        blanks - str
            String with n blanks
    """
    return ' '*n

def dict_products(dicts):
    """
    Gives every combination over a dictionary that has values of some iteratable object

    Parameters
    ----------
        dicts - dict
            Dictionary that has values that are iteratable.
            e.g. dicts = {
                    'empty': [None]
                    'even_number': [2, 4],
                    'letter': ('a', 'b', 'c'),
                    }
    Returns
        dict_combinations - tuple of dicts
            For the example above, will return:
                (
                    {'empty': None, 'even_number': 2, 'letter': 'a'},
                    {'empty': None, 'even_number': 2, 'letter': 'b'},
                    {'empty': None, 'even_number': 2, 'letter': 'c'},
                    {'empty': None, 'even_number': 4, 'letter': 'a'},
                    {'empty': None, 'even_number': 4, 'letter': 'b'},
                    {'empty': None, 'even_number': 4, 'letter': 'c'},
                )
    """
    return (dict(zip(dicts, x)) for x in itertools.product(*dicts.values()))

def interpolate(x_low, x_high, y_low, y_high, x):
    """
    Linear interpolation

    Parameters
    ----------
        x_low - float
            Lower x value for interpolation
        x_high - float
            Higher x value for interpolation
        y_low - float
            Lower y value that corresponds to x_low for interpolation
        y_high - float
            Higher y value that corresponds to x_high for interpolation
        x - float
            Target x value
    Returns
    -------
        y - float
            Target y value
    """
    if x == x_low:
        return y_low
    elif x == x_high:
        return y_high

    m = (y_high - y_low)/(x_high - x_low)
    b = y_high - m * x_high
    return m * x + b

def saveplt(ax, filename):
    """
    Save matplotlib figure as pickle file

    Parameters
    ----------
        ax - matplotlib.axes.Axes object
            Matplotlib object to save
        filename - string
            Name of file to save the figure to
    """
    with open(filename, 'wb') as f_ptr:
        pickle.dump(ax, f_ptr)

def loadplt(filename, show = False):
    """
    Save matplotlib figure as pickle file

    Parameters
    ----------
        filename - string
            Name of file to save the figure to
        show - boolean
            Whether or not to show the plot
    Returns
    -------
        ax - matplotlib.axes.Axes object
    """
    from matplotlib import pyplot as plt

    with open(filename, 'rb') as f_ptr:
        ax = pickle.load(f_ptr)
    if show:
        plt.show()
    else:
        return ax