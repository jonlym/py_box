# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 07:52:17 2016

@author: Jonathan Lym
"""

import numpy as np


def fit_CpoR(T, CpoR):
    mid = 800.
    mid_index = np.where(T >= mid)[0][0]
    T_low = T[:mid_index]
    CpoR_low = CpoR[:mid_index]
    
    T_high = T[mid_index:]    
    CpoR_high = CpoR[mid_index:]
    
    a_low_rev = np.polyfit(x = T_low, y = CpoR_low, deg = 4)
    a_high_rev = np.polyfit(x = T_high, y = CpoR_high, deg = 4)

    
    empty_arr = np.array([0]*2)
    a_low = np.concatenate((a_low_rev[::-1], empty_arr))
    a_high = np.concatenate((a_high_rev[::-1], empty_arr))
    return (a_low, a_high, T[mid_index])