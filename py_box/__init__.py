# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 21:10:42 2016

@author: Jon Lym
"""

def any_alpha(string):
    """Returns True if any alphabetic characters are in the string. False otherwise"""
    for character in string:
        if character.isalpha():
            return True
    else:
        return False
