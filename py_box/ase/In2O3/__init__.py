# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 17:12:09 2017

@author: Jonathan Lym
"""

from ase.io import read, write
import heapq
import numpy as np
import os

configuration_template = """import numpy as np
from ase.io import read
from ase.io.trajectory import Trajectory
from ase.optimize import LBFGS
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms
from os.path import basename, exists
from py_box.ase import run_testRun
from py_box.ase.In2O3 import index_dict
from py_box.ase.set_calc import print_vasp_param, calc_dict, handle_restart
from py_box.cluster_expansion import run_cluster_expansion

testRun = False
file_name_py = basename(__file__)
file_name = file_name_py.replace('.py','')
file_traj = file_name + '.traj'
file_out = file_name +'.out'
start_file = './POSCAR_start'

#Default calculator
vasp_param = calc_dict['In2O3']
#Modifications to calculator
vasp_param['setups'] = {'In': '_d'}
vasp_param['nelm'] = 400
vasp_param['nelmdl'] = -10
vasp_param['lwave'] = False

try:
    #Check for a restarted calculation
    sys = Trajectory(file_traj, 'r')[-1]
except (IOError, RuntimeError):
    #Read from clean surface
    print "Importing trajectory file from: %s" % start_file
    sys = read(start_file)
    mode = 'w'
else:
    print "Importing trajectory file from: %s" % file_traj
    mode = 'a'

calc = Vasp(**vasp_param)
#Checking for geometric optimization or geometric optimization
if exists(file_out):
    if 'Starting energy optimization.' not in open(file_out, 'r').read():
        calc.set(kpts = (1, 1, 1), ispin = 1)
else:
    calc.set(kpts = (1, 1, 1), ispin = 1)
handle_restart(calc, sys)
print_vasp_param(calc)
sys.set_calculator(calc)

del sys.constraints
print 'Constraints:'
print '	c1: Fixed lowest layer'
c1 = FixAtoms(mask = [a.z < 12.1 for a in sys])
sys.set_constraint(c1)

if testRun:
    run_testRun(sys)
else:
    geo_traj = Trajectory(file_traj, mode, sys)
    dyn = LBFGS(sys, trajectory = geo_traj)
    if calc.input_params['kpts'] == (1, 1, 1):
        print 'Starting geometric optimization.'
        dyn.run(fmax = 0.05)
        print 'Completed geometric optimization.'
        calc = Vasp(**vasp_param)
        handle_restart(calc, sys)
        sys.set_calculator(calc)
        print_vasp_param(calc)
    print 'Starting energy optimization.'
    dyn.run(fmax = 0.05)
    print 'Completed energy optimization'
    energy = sys.get_potential_energy()
    print "Energy: %f" % energy
print "Finding next configuration..."
run_cluster_expansion(train_path = '/project/projectdirs/m1893/jlym/In2O3/unit_cell_relaxed/cluster_expansion/train', clusters_path = '/project/projectdirs/m1893/jlym/In2O3/unit_cell_relaxed/cluster_expansion/cluster_vacancy.xlsx', configs_all_path = '/project/projectdirs/m1893/jlym/In2O3/unit_cell_relaxed/cluster_expansion/config_vacancy.xlsx', log_path = '/project/projectdirs/m1893/jlym/In2O3/unit_cell_relaxed/cluster_expansion/output.log', submit_job = testRun)
print "Completed %s" % file_name_py"""

index_dict = {'In1': [26, 27],
              'In2': [24, 25],
              'In3': [15, 14],
              'In4': [22, 23],
              'O1':  [75, 74],
              'O2':  [76, 77],
              'O3':  [46, 47],
              'O4':  [73, 72],
              'O5':  [43, 42],
              'O6':  [68, 69],
              'In1A': 26,
              'In1B': 27,
              'In2A': 24,
              'In2B': 25,
              'In3A': 15,
              'In3B': 14,
              'In4A': 22,
              'In4B': 23,
              'O1A':  75,
              'O1B':  74,
              'O2A':  76,
              'O2B':  77,
              'O3A':  46,
              'O3B':  47,
              'O4A':  73,
              'O4B':  72,
              'O5A':  43,
              'O5B':  42,
              'O6A':  68,
              'O6B':  69}

def find_closest_H(atoms, i, n = 1):
    """Finds the 'n' indices of the H atoms closest to atom 'i'."""
    all_H_indices = []
    all_H_dist = []
    for j, atom in enumerate(atoms):
        if atom.symbol == 'H':
            all_H_indices.append(j)
            all_H_dist.append(atoms.get_distance(i, j))
    H_low_dist = heapq.nsmallest(n, all_H_dist)
    H_low_indices = []
    for dist in H_low_dist:
        H_low_indices.append(all_H_indices[all_H_dist.index(dist)])
    return H_low_indices

def get_new_index(atoms_indices, vacancy_indices, start_index = 0):
    """
    Returns the adjusted indices given the vacancies already created on the surface.
    start_index = 1 should be used when the indices are not 0 indexed.
    """
    out_atom_indices = []

    if vacancy_indices is None:
        return atoms_indices
    elif type(vacancy_indices) is int:
        vacancy_indices = [vacancy_indices]

    if type(atoms_indices) is int:
        atom_indices = [atoms_indices]

    for atom_index in atom_indices:
        if atom_index in vacancy_indices:
            out_atom_indices.append(np.nan)
        else:
            offset = start_index - sum([1 for vacancy_index in vacancy_indices if vacancy_index < atom_index])
            out_atom_indices.append(atom_index + offset)

    if len(out_atom_indices) == 1:
        return out_atom_indices[0]
    else:
        return out_atom_indices

def get_In2O3_configuration(bitstring, width = 12):
    indices = [75, 76, 46, 73, 43, 68, 74, 77, 47, 72, 42, 69]
    del_indices = []
    for i, j in enumerate(bitstring):
        if j == '0':
            del_indices.append(indices[i])
    cluster = os.environ['CLUSTER']
    if 'farber' in cluster:
        clean_path = '/home/work/ccei_biomass/users/jlym/In2O3/unit_cell_relaxed/In2O3_110_clean/CONTCAR'
    elif 'edison' in cluster:
        clean_path = '/project/projectdirs/m1893/jlym/In2O3/unit_cell_relaxed/In2O3_110_clean/CONTCAR'
    elif 'work' in cluster:
        clean_path = 'C:\\Users\\Jonathan Lym\\Google Drive\\UDel Documents\\UDel Research\\In2O3\\D Orbitals\\In2O3_110_clean\\CONTCAR'
    else:
        raise Exception("Path to In2O3 clean not specified.")
    atoms = read(clean_path)
    del atoms[del_indices]
    return atoms

def run_In2O3_configuration(configuration, job_array = False, job_file = 'joblist.txt', folder_file = 'folderlist.txt', rel_path = None, submit_job = True):
    """
    Writes the In2O3 ASE script and adds to the Job Array files
    :param configuration:
    :param job_file:
    :param folder_file:
    :param rel_path:
    :return: A boolean indicating whether or not a job has been submitted
    """

    home_dir = os.getcwd()
    #If path not specified, make a folder in the current directory with that name
    if rel_path is None:
        rel_path = './{}/'.format(configuration.name)

    #If the folder does not exist, make it
    if not os.path.isdir(os.path.join(rel_path, configuration.name)):
        os.makedirs(os.path.join(rel_path, configuration.name))
    elif os.path.isfile(os.path.join(rel_path, '{}/POSCAR_start'.format(configuration.name))):
        return False
    #Write the configuration to a file
    write(os.path.join(rel_path, '{}/POSCAR_start'.format(configuration.name)), get_In2O3_configuration(bitstring = convert_sigma_to_bitstring(configuration.sigma)))

    #Write the script
    script_name = '{}.py'.format(configuration.name)
    with open(os.path.join(rel_path, '{}/{}'.format(configuration.name, script_name)), 'w') as py_ptr:
        py_ptr.write(configuration_template)

    if job_array:
        #Add the file to the job array
        with open(job_file, 'a') as job_ptr:
            job_ptr.write(configuration.name)
        with open(folder_file, 'a') as folder_ptr:
            folder_ptr.write('{}{}'.format(rel_path, configuration.name))
    else:
        os.chdir(os.path.join(rel_path, configuration.name))
        if submit_job:
            cluster = os.environ['CLUSTER']
            if 'farber' in cluster:
                os.system('qase vasp 20 {} -s'.format(script_name))
            elif 'edison' in cluster:
                os.system('qvasp_ase_log 24 {} -w 10:00:00 -s'.format(script_name))
        os.chdir(home_dir)
    return True

def convert_sigma_to_bitstring(sigma):
    """
    Converts a sigma vector to a bitstring
    [1, -1, 1, 1] becomes '1011'
    :param sigma: Configuration vector containing -1 and 1 separated by commas
    :return:
    """
    if type(sigma) is np.ndarray or type(sigma) is list:
        bitstring = ''
        for i in sigma:
            if i == -1:
                bitstring = '{}0'.format(bitstring)
            else:
                bitstring = '{}1'.format(bitstring)
        return bitstring
    elif type(sigma) is str:
        return sigma.replace('-1', '0').replace(',', '').replace(' ', '')
    else:
        print type(sigma)
        raise Exception('Invalid data type when converting sigma to bitstring.')
