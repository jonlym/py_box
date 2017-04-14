# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 16:02:35 2017

@author: Jonathan Lym
"""

import numpy as np
import matplotlib.pyplot as plt

def plot_dos(dos_file = 'dos_output.txt', spin = True):
    if spin:
        orbitals = ['s up', 's down', 'p up', 'p down', 'd up', 'd down', 'f up', 'f down']
        line_types = ['r-', 'r--', 'b-', 'b--', 'k-', 'k--', 'g-', 'g--']
    else:
        orbitals = ['s', 'p', 'd', 'f']
        line_types = ['r-', 'b-', 'k-', 'g-']
    with open(dos_file, 'r') as dos_ptr:
        lines = dos_ptr.readlines()

        #Determine the number of columns
        buf = lines[0].split(" ")
        buf = [element for element in buf if element != '']
        n_col = len(buf)-1
        
        energies = np.zeros(shape = len(lines))
        dos = np.zeros(shape = (n_col, len(lines)))
        #Extract data
        for i, line in enumerate(lines):
            line = line.replace('\n', '')
            buf = line.split(" ")
            buf = [element for element in buf if element != '']
            energies[i] = float(buf[0])
            for j, buf_data in enumerate(buf[1:]):
                if 'down' in orbitals[j]: 
                    dos[j, i] = -float(buf_data)
                else:
                    dos[j, i] = float(buf_data)
        plt.figure()
        for i, dos_orbital in enumerate(dos):
            plt.plot(energies, dos_orbital, line_types[i])
            plt.legend(orbitals[:n_col])
            plt.ylabel('Density of States')
            plt.xlabel('Energy (eV)')
        
        ax = plt.gca()
        ylim = ax.get_ylim()
        plt.xlim([min(energies), max(energies)])
        plt.plot([0, 0], ax.get_ylim(), '#808080')
        plt.ylim(ylim)        