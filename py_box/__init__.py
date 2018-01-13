# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 21:10:42 2016

@author: Jon Lym
"""

import numpy as np
from itertools import product, izip
from datetime import datetime


def base10_to_base3(n, width):
    nums = np.zeros(shape = (width, ))
    if n == 0:
        return nums
    i = 0
    while n:
        n, r = divmod(n, 3)
        nums[i] = r
        i += 1
    return nums[::-1]

def any_alpha(string):
    """Returns True if any alphabetic characters are in the string. False otherwise"""
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
    return keys.keys()

def plot_parity(x, y, decimals = 2):
    from matplotlib import pyplot as plt
    data = np.concatenate((x, y), axis = 0)
    min_data = np.round(min(data*10.**decimals))/10.**decimals
    max_data = np.round(max(data*10.**decimals))/10.**decimals
    fig = plt.figure()
    plt.plot(np.array([min_data, max_data]), np.array([min_data, max_data]), 'k-', x, y, 'bo')
    axes = plt.gca()
    axes.set_ylim([min_data, max_data])
    axes.set_xlim([min_data, max_data])
#    plt.plot(np.array([min_data, min_data]), np.array([max_data, max_data]), 'k-')
    return (fig, axes)

def get_time():
    return str(datetime.now())

def get_null(mat, rtol=1e-5):
    u, s, v = np.linalg.svd(mat)
    rank = (s > rtol*s[0]).sum()
    return v[rank:].T.copy()

def get_RMSE(xs_data, xs_fit):
    return np.sqrt(np.mean([(x_data-x_fit)**2 for x_data, x_fit in zip(xs_data, xs_fit)]))

def spherical_to_xyz(r = 1., theta = 0., phi = 0., degrees = True):
    if degrees:
        theta = np.radians(theta)
        phi = np.radians(phi)
    x = r * np.cos(theta) * np.sin(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(phi)
    return np.array([x, y, z])

def get_n_blanks(n):
    return ' '*n

def dict_product(dicts):
    return (dict(izip(dicts, x)) for x in product(*dicts.itervalues()))