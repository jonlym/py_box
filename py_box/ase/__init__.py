# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 16:38:29 2016

@author: Jonathan Lym
"""

from platform import system
from pprint import pprint
from ase.visualize import view

def run_testRun(atom_obj):
    print 'Test Run'
    os_name = system()
    if type(atom_obj) is list:        
        for i in range(len(atom_obj)):    
            pprint(vars(atom_obj[i]))
    else:
        pprint(vars(atom_obj))
    if os_name.lower() == 'linux':
        view(atom_obj)