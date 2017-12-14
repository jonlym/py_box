# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 21:10:42 2016

@author: Jon Lym
"""

import numpy as np
from datetime import datetime


def base10_to_base3(n: int, width: int):
    """
    Converts base 10 numbers to base 3 width.
    """
    nums = np.zeros(shape = (width, ))
    if n == 0:
        return nums
    i = 0
    while n:
        n, r = divmod(n, 3)
        nums[i] = r
        i += 1
    return nums[::-1]

def any_alpha(string: str):
    """
    Returns True if any alphabetic characters are in the string. False otherwise.
    """
    for character in string:
        if character.isalpha():
            return True
    else:
        return False

def get_unique_list(data: list):
    """
    Given a list, returns a unique list.
    """
    keys = {}
    for item in data:
        keys[item] = 1
    return keys.keys()

def plot_parity(x: np.ndarray, y: np.ndarray, decimals: int = 2):
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

def get_null(mat: np.ndarray, rtol: float=1.e-5):
    """
    Returns the nullspace of a 2D matrix, mat
    """
    u, s, v = np.linalg.svd(mat)
    rank = (s > rtol*s[0]).sum()
    return v[rank:].T.copy()

def get_RMSE(xs_data: np.ndarray, xs_fit: np.ndarray):
    """
    Returns the root mean squared error given two vectors of equal lengths
    """
    return np.sqrt(np.mean([(x_data-x_fit)**2 for x_data, x_fit in zip(xs_data, xs_fit)]))

def spherical_to_xyz(r: float = 1., theta: float = 0., phi: float = 0., degrees: bool = True):
    """
    Converts spherical coordinates to Cartesian coordinates. Angles are in degrees by default. Set degrees to False
    to use radians
    """
    if degrees:
        theta = np.radians(theta)
        psi = np.radians(phi)
    x = r * np.cos(theta) * np.sin(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(phi)
    return np.array([x, y, z])

def get_n_blanks(n: int):
    return ' '*n
