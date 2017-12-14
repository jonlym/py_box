# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 12:25:37 2017

@author: Jonathan Lym
"""

from datetime.datetime import now
from os.path import exists

def update_log(file_name, txt):
    log_file = file_name + '.log'
    if exists(log_file):
        mode = 'a'
    else:
        mode = 'w'
        
    with open(log_file, mode) as log_ptr:
        log_ptr.write('%s  %s\n' (now(), txt))