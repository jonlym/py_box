# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 08:54:32 2016

@author: Jonathan Lym
"""

import numpy as np
import shutil as sh
import re
from ase.io import read, write
from ase.visualize import view
from ase.neb import interpolate
from os import chdir, getcwd, system, makedirs, walk
from os.path import relpath, expanduser, exists, join
from py_box.ase.TS import label_folder, compare_initial_final

def submit_CM(n_images = 6, nest = 1, wall_time = '12:00:00', version = 'v54', queue = None, submit_jobs = True):
    """
    Submits the constrained minimization python files to the cluster using the qase script developed by Glen Jeness.
    Parameters supported:
        n_images: Number of intermittent images in folders ## where # is a number. Default is 6.
        nest: Indicates how nested the trial is in the constrained minimization folder. Default is 1.
        wall_time: Wall time on Edison in format: HH:MM:SS. Default is 12:00:00 (12 Hours).
        version: Compatible with Farber and Edison. Default is v54
        queue: Compatible with Edison. Enter 'debug' to enter the debug queue. Default is None
        submit_jobs: Can be toggled to automatically submit qs files or not. Default is True.
    In order to detect the cluster, 'misc_info.txt' should be in the home directory and should contain the name of the cluster (e.g. squidward, farber, edison)
    """
    file_base = relpath('.', '../'*nest).replace('/', '_')
    email = 'jlym@udel.edu'
    print '-'*20
    print "Files to be submitted take the form: %s##.py where ## range between 01 and %s" % (file_base, label_folder(n_images))
    if not submit_jobs:
        print "WARNING: submit_jobs set to False. Jobs will not be submitted!"
        
    #Determines the cluster to use the appropriate submit command
    home = expanduser('~')
    cluster_file = open('%s/misc_info.txt' % home, 'r')
    content = cluster_file.read()
    cluster_file.close()
    if 'squidward' in content.lower():
        cluster = 'squidward'
        n_cores = 16
        qase_start = 'qase' 
        qase_options = ''
        if submit_jobs:
            qase_options = '%s -s' % qase_options
    
    elif 'farber' in content.lower():
        cluster = 'farber'
        n_cores = 20
        qase_start  = 'qase vasp'
        if submit_jobs:
            qase_options = '%s -s' % (qase_options, version)
        if queue is not None:
            qase_options = '%s -q %s' % (qase_options, queue)
        qase_options = '-p %s' % version
        
    elif 'edison' in content.lower():
        cluster = 'edison'
        qase_start = 'qvasp_ase'
        qase_options = '%s -w %s -p %s' % (qase_options, wall_time, version)
        if submit_jobs:
            qase_options = '%s -s' % qase_options
        if queue is not None:
            if queue is 'knl':
                qase_options = '%s -q %s' % (qase_options, queue)
                n_cores = 32
                cluster = 'cori'
            else:
                qase_options = '%s -q %s' % (qase_options, queue)
                n_cores = 24

    else:
        print 'Warning. None of the compatible cluster types found in misc_info.txt' 

    #Information related to job
    print 'Cluster: %s' % cluster
    print 'Cores per job: %d' % n_cores
    if 'farber' in cluster or 'edison' in cluster:
        print 'Using Vasp version %s' % version 
        if 'edison' in cluster or 'cori' in cluster:
            print 'Walltime per job: %s' % wall_time
    print '-'*20

    for i in range(1, n_images+1):
        print "Processing %i" % i
        folder = label_folder(i)
        chdir(folder)
        system('%s %d %s%s.py %s' % (qase_start, n_cores, file_base, folder, qase_options))        
        
        if 'squidward' in cluster:
            print 'Adding e-mail notification to file'
            qs_file = open('%s%s.qs' % (file_base, folder), 'r')
            lines = qs_file.readlines()
            qs_file.close()
            for i, line in enumerate(lines):
                if '#$' not in line and '#!' not in line: #At the end of bash options
                    lines.insert(i, '#$ -m beas\n#$ -M %s\n' % email)
                    break
            qs_file = open('%s%s.qs' % (file_base, folder), 'w')
            lines = "".join(lines)
            qs_file.write(lines)
            qs_file.close()

            if submit_jobs:           
                system('qsub %s%s.qs' % (file_base, folder))
        chdir('..')
    print "Completed submit_NEB"

def initialize_CM(n_images = 8, nest = 1, write_POSCAR = True, write_python = True, python_template = 'template.py'):
    """
    Initializes the NEB calculation.
    Parameters supported:
        n_images: Number of intermittent images in folders ## where # is a number. Default is 6.
        nest: Indicates how nested the trial is in the constrained minimization folder. Default is 1.
        write_POSCAR: If True, it will write a POSCAR file to the folders as POSCAR_start. Default is True.
        write_python: If True, it will copy the template file to the folders. Default is True.
    """
        
    file_base = relpath('.', '../'*nest).replace('/', '_')
    
    #Read initial and final state
    initial = read('./initial/CONTCAR')
    final = read('./final/CONTCAR')
    
    #Generating the interpolated images
    images = [initial]
    for i in range(n_images):
        images.append(initial.copy())
    images.append(final)
    interpolate(images)
    
    #for i in range(1, n_images+1):
    #    print 'Checking bond lengths for image %d' % i
    #    check_bond_lengths(images[i], False)
    
    compare_initial_final(initial, final)
    
    for i in range(0, n_images+2):
        print 'Processing image %d' % i
        folder = label_folder(i)
        if not exists(folder):
            makedirs(folder)
        if write_POSCAR:
            POSCAR_path = join(folder, 'POSCAR_start')
            print 'Writing image %d to %s' % (i, POSCAR_path)
            write(POSCAR_path, images[i])
        if write_python:
            python_path = join(folder, '{}{}.py'.format(file_base, folder))
            print 'Copying template.py to %s' % python_path
            sh.copyfile(python_template, python_path)
    print 'Completed initialize_NEB'

