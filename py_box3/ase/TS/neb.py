# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 08:54:32 2016

@author: Jonathan Lym
"""

import numpy as np
import re
from ase.io import read, write
from ase.visualize import view
from ase.neb import interpolate as ase_interpolate
from shutil import copy, copytree
from os import chdir, getcwd, system, makedirs, walk
from os.path import relpath, expanduser, exists, join, isfile

def compare_initial_final(initial_image, final_image):
    """
    Compares the relative locations of the initial and the final image. Also
    checks if the indices of the two images represent the same elements.
    """
    initial_len = len(initial_image)
    final_len = len(final_image)
    if initial_len == final_len:
        print("Initial image and final image contain the same number of atoms.")
        min_len = initial_len
    elif initial_len > final_len:
        print("Initial image has more atoms than final image.")
        min_len = final_len
    else:
        print("Final image has more atoms than initial image.")
        min_len = initial_len

    print("Index\tInitial\tFinal\tDelta r (A)\tDelta x (A)\tDelta y (A)\tDelta z (A)")
    print("-----------------------------------------------------------------")
    for i in range(min_len):
        del_x = final_image[i].x - initial_image[i].x
        del_y = final_image[i].y - initial_image[i].y
        del_z = final_image[i].z - initial_image[i].z
        del_r = np.sqrt(del_x**2 + del_y**2 + del_z**2)
        print(("%d\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f" % (i, initial_image[i].symbol, final_image[i].symbol, del_r, del_x, del_y, del_z)))
        if initial_image[i].symbol != final_image[i].symbol:
            print("Symbol mismatch!")
    print("-----------------------------------------------------------------")

def move_atoms_by_cell(in_path, atoms_list, cell_dir_list, offset_list):
    """
    Opens the file based on the file path and moves the atoms marked in the 
    atoms_list with the indices provided by cell_dir by a multiple offset.
    
    Positions will be moved as:
    New Position = Old Position + offset * cell[cell_dir]
    """
    atoms_obj = read(in_path)
    for (i, cell_i, offset) in zip(atoms_list, cell_dir_list, offset_list):
        cell_dir = offset*atoms_obj.get_cell()[abs(cell_i)]
        print(("Moving atom %s %d by %s" % (atoms_obj[i].symbol, i, np.array_str(cell_dir))))
        atoms_obj[i].position += cell_dir
    view(atoms_obj)
    print(("Overwrite %s (Y/N)?" % in_path))
    overwrite = eval(input('>>'))
    if overwrite.lower() == 'y':
        out_path = in_path
    else:
        out_path = in_path + '_fix'
    write(out_path, atoms_obj)
    print("Completed ase_NEB.move_atoms_by_cell")
    print(("Result written to %s" % out_path))

def label_folder(num):
    """
    Takes the number and returns the string name in the format '##'
    """
    if num > 9:
        return str(num)
    else:
        return '0'+str(num)

def initialize_NEB(location = '.'):
    """
    Uses the template folder (expected to be in one directory up) to initialize NEB calculation.
    """
    chdir(location)
    dir_name = relpath('.', '..')
    if exists('../template'):
        print('Template folder found.')
        if isfile('../template/python_template.py'):
            print(('Python script found. Copying to location as %s.py' % dir_name))
            copy('../template/python_template.py', './%s.py' % dir_name)
        if exists('../template/initial'):
            print('Initial folder found. Copying its contents.')
            copytree('../template/initial', './initial')
        if exists('../template/final'):
            print('Final folder found. Copying its contents.')
            copytree('../template/final', './final')

def print_energies():
    """
    Tabulates the energies from the OUTCAR files. If print_rel_energy is set to True, the initial energy will be found in initial/OUTCAR and the table will
    have an extra field showing energies relative to this value.
    """
    #Setting up np.array type
    values = []
    dtype = [('path', 'a40'), ('energy', float)]
    print("nan indicates the run is not finished.")
    print("inf indicates the run does not exist.")
    #Try reading the initial image for relative energies
    try:
        initial = read('initial/OUTCAR')
    except (IndexError, IOError):
        print("Initial image not found.")
        print_rel_energy = False
    else:
        initial_energy = initial.get_potential_energy()
        print_rel_energy = True

    #Print table headings:
    if print_rel_energy:
        print("Path\t\tEnergy(eV)\tEnergy relative to image 0 (eV)")
        print(("-"*30))
        print(("initial\t\t%.5f\t%.2f" % (initial_energy, initial_energy-initial_energy)))
    else:
        print("Path\t\tEnergy(eV)")
        print(("-"*30))

    for root, dirs, files in walk(getcwd()):
        for directory in dirs:
            #If the folder has the form ## but not 00
            if re.search('([0-9][1-9]|[1-9][0-9])', directory) is not None and '-' not in directory:
                #Check if folder is interpolated.
                full_path = join(root, directory)
                image_path = re.search('[0-9]{2}-[0-9]{2}.*([0-9][1-9]|[1-9][0-9])', full_path)
                if image_path is None:
                    image_path = directory
                else:
                    image_path = image_path.group(0)

                #Find energy
                try:
                    image = read('%s/OUTCAR' % full_path)
                except IndexError:
                    energy = np.nan
                except IOError:
                    energy = np.inf
                else:
                    energy = image.get_potential_energy()
                values.append((image_path, energy))
    #Sort values and print
    unsorted_values = np.array(values, dtype=dtype)
    sorted_values = np.sort(unsorted_values, order = 'path')
    for values in sorted_values:
        if print_rel_energy:
            print(("%s\t\t%.5f\t%.2f" % (values['path'], values['energy'], values['energy'] - initial_energy)))
        else:
            print(("%s\t\t%.5f" % (values['path'], values['energy'])))

    try:
        final = read('final/OUTCAR')
    except (IndexError, IOError):
        print("Final image not found")
    else:
        if print_rel_energy:
            final_energy = final.get_potential_energy()
            print(("final\t\t%.5f\t%.2f" % (final_energy, final_energy-initial_energy)))
        else:
            final_energy = final.get_potential_energy()
            print(("final\t\t%.5f" % final_energy))

def interpolate(initial_atoms, final_atoms, n_images):
	"""
	Interpolates the positions between initial_atoms and final_atoms

	Parameters
	----------
		initial_atoms : ase.Atoms
			Initial position of atoms
		final_atoms : ase.Atoms
			Final position of atoms
		n_images : int
			Number of images (not including initial_atoms and final_atoms)
	Returns
	-------
		images : list of ase.Atoms
			Interpolated images (the 0th index is the initial image and the -1th index is the final image)
	"""
	initial_copy = initial_atoms.copy()
	images = [initial_copy]
	for i in range(n_images):
		images.append(initial_copy.copy())
	images.append(final_atoms)
	ase_interpolate(images)
	return images