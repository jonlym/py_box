# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 21:10:42 2016

@author: Jon Lym
"""

import itertools
import numpy as np
from datetime import datetime
import pickle
from matplotlib import pyplot as plt

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
    with open(filename, 'rb') as f_ptr:
        ax = pickle.load(f_ptr)
    if show:
        plt.show()
    else:
        return ax