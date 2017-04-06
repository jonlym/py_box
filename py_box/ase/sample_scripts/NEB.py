# -*- coding: utf-8 -*-
"""
Created on Thu Nov 03 12:31:25 2016

@author: Jonathan Lym
"""

from ase.io import read
from ase.io.vasp import write_vasp
from ase.io.trajectory import Trajectory as traj
from ase.neb import interpolate
from ase.constraints import FixAtoms
import shutil as sh
import os
from os.path import basename
from run_testRun import run_testRun
from set_calc import set_neb_calc_In2O3
from ase_TS import label_folder


testRun = False
restart = False

#NEB Parameters
if restart:
    #Use high convergence
    encut = 400
    ispin = 1
    kpts = (4, 3, 1)
    ediffg = -0.05
else:
    #Use low convergence
    encut = 400
    ispin = 1
    kpts = (4, 3, 1)
    ediffg = -0.05
NIMAGES = 9


print "Restart = %r" % restart
print "NIMAGES = %d" % NIMAGES

#Read initial and final state
initial = read('../In2O3_110_H_24/In2O3_110_H_24.traj')
final_path = '../In2O3_110_H_76/In2O3_110_H_76.traj'
final = read(final_path)

#Generating the interpolated images
images = [initial]
for i in range(NIMAGES):
    images.append(initial.copy())
images.append(final)
interpolate(images)

for i in range(0, (NIMAGES+2)):
    print 'Processing image %d' % i
    dir = label_folder(i)
    if not os.path.exists(dir):
        os.makedirs(dir)
    if restart:
        if i != 0 and i != (NIMAGES+1):
            print "Restarting calculation. Copying CONTCAR to POSCAR"
            image = read(dir+'/CONTCAR')
        else:
            image = read(dir+'/POSCAR')
    else:
        print "Starting calculation. Writing POSCAR from interpolation"
        image = images[i]
        #Place adjustments here
    calc = set_neb_calc_In2O3(atoms_obj = image, NIMAGES = NIMAGES, kpts = kpts, encut = encut, ispin = ispin, ediffg = ediffg)
    print "Constraints:"
    print "\tc1: Freeze 2 lower layers"
    c1 = FixAtoms(mask = [atom.z < 8.5 for atom in image])
    image.set_constraint(c1)
    if not testRun:
        calc.initialize(image)
        write_vasp(dir+'/POSCAR', calc.atoms_sorted, symbol_count=calc.symbol_count)
        calc.write_potcar()
print "Successfully written files"
calc = set_neb_calc_In2O3(initial, NIMAGES, kpts = kpts, encut = encut, ispin = ispin, ediffg = ediffg)
initial.set_calculator(calc)
if testRun:
    run_testRun(images)
else:
    print "Running NEB"
    print initial.get_potential_energy()
print "Completed %s" % basename(__file__)
#    if restart:
#	print "Constraints:"
#        print "\tc1: Fixing all surface atoms"
#        for i in range(1, NIMAGES+1):
#            print 'Processing image %d' % i
#            dir = label_folder(i)
#            if os.path.exists(dir) == False:
#                os.makedirs(dir)
#            sh.copyfile(dir+'/CONTCAR', dir+'/POSCAR')
#            atoms = read(dir+'POSCAR')
#            c1 = FixAtoms(mask=[atom.symbol != 'H' for atom in atoms])
#            atom.
#    else:
#        for i in range(0, (NIMAGES+2)):
#            print 'Processing image %d' % i
#            dir = label_folder(i)
#            if os.path.exists(dir) == False:
#                os.makedirs(dir)
#            
#            image = images[i]
#            calc = set_neb_calc(image, NIMAGES)
#            calc.initialize(image)
#            
#            write_vasp(dir+'/POSCAR', calc.atoms_sorted, symbol_count=calc.symbol_count)
#            calc.write_potcar()
#        print "Successfully written POSCAR files"
#    calc = set_neb_calc(initial, NIMAGES)           
#    print "Calculating energies..."
#    initial.set_calculator(calc)
#    print initial.get_potential_energy()
#print 'Completed %s' % basename(__file__)
