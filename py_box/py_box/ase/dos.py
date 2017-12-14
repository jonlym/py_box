# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 16:02:35 2017

@author: Jonathan Lym
"""

import warnings
import numpy as np
import collections

column_dict = {'Energy': 0,
               's': 1,
               'p': 2,
               'd': 3,
               'f': 4,
               's up': 1,
               's down': 2,
               'p up': 3,
               'p down': 4,
               'd up': 5,
               'd down': 6,
               'f up': 7,
               'f down': 8}


def read_dos_output(dos_file = 'dos_output.txt', spin = True, orbitals = ['s up', 's down', 'p up', 'p down', 'd up', 'd down', 'f up', 'f down']):
    """Read the dos output file in the format used by Dr. Glen Jenness' script."""
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
    dos_data = collections.namedtuple('dos_data', ['energies', 'dos'])
    return dos_data(energies = energies, dos = dos)

def plot_dos(dos_file = 'dos_output.txt', spin = True):
    """Plots the density of states given a dos output file in the format used by Dr. Glen Jenness' script."""
    import matplotlib.pyplot as plt

    if spin:
        orbitals = ['s up', 's down', 'p up', 'p down', 'd up', 'd down', 'f up', 'f down']
        line_types = ['r-', 'r--', 'b-', 'b--', 'k-', 'k--', 'g-', 'g--']
    else:
        orbitals = ['s', 'p', 'd', 'f']
        line_types = ['r-', 'b-', 'k-', 'g-']
    dos_data = read_dos_output(dos_file = dos_file, spin = spin, orbitals = orbitals)
    plt.figure()
    for i, dos_orbital in enumerate(dos_data.dos):
        plt.plot(dos_data.energies, dos_orbital, line_types[i])
    plt.legend(orbitals[:np.shape(dos_data.dos)[0]])
    plt.ylabel('Density of States')
    plt.xlabel('Energy (eV)')
    
    ax = plt.gca()
    ylim = ax.get_ylim()
    plt.xlim([min(dos_data.energies), max(dos_data.energies)])
    plt.plot([0, 0], ax.get_ylim(), '#808080')
    plt.ylim(ylim)        
        
def get_mean_energy(dos_file = 'dos_output.txt', spin = True, orbital = 'd up', interval = [-30, 30]):
    dos_data = read_dos_output(dos_file = dos_file, spin = spin)
    i_min = get_nearest(dos_data.energies, interval[0])
    i_max = get_nearest(dos_data.energies, interval[1])
#    print 'Interval: {} to {}'.format(interval[0], interval[1])
    energies = dos_data.energies[i_min:(i_max+1)]
    dos_data_orbital = dos_data.dos[column_dict[orbital]-1][i_min:(i_max+1)]
    weighted_data = energies * dos_data_orbital
#    print 'Numerator: {}'.format(integrate(energies, weighted_data))
#    print 'Denominator: {}'.format(integrate(energies, dos_data_orbital))
    m_energy =  integrate(energies, weighted_data)/integrate(energies, dos_data_orbital)
    return m_energy

def integrate(x, y):
    """Numerically integrates the data using the trapesium method."""
    area = 0.
    for i in range((len(x)-1)):
        area += 0.5 * (x[i+1] - x[i]) * (y[i+1] + y[i])
    return area

def get_nearest(array, value):
    """Finds the index closest to the value."""
    return (np.abs(array-value)).argmin()
