# -*- coding: utf-8 -*-
"""
Created on Fri Mar 03 14:36:24 2017

@author: Jon Lym
"""

from ase.io import read
from ase.constraints import FixAtoms
from set_calc import set_dimer_calc_ZnOCu
from run_testRun import run_testRun
import numpy as np
import shutil
import glob

def find_last_modecar():
    max_val = -1
    print((glob.glob('./MODECAR*')))
    for modecar in glob.glob('./MODECAR*'):
        try:
            i = int(float(modecar[9:]))
        except ValueError:
            continue
        if i > max_val:
            max_val = i
    return max_val

restart = True
testRun = False
TS_path = '../03/OUTCAR'
TS_vector_path = '../02/OUTCAR'

if not restart:
    TS = read(TS_path)
    TS_vector = read(TS_vector_path)
    print("Creating MODECAR using:")
    print(("%s as TS" % TS_path))
    print(("%s as approximation for MODECAR" % TS_vector_path))
    modecar = np.array(TS_vector.get_positions() - TS.get_positions())
    np.savetxt('MODECAR', modecar, fmt = '%17.8f')
else:
    print("Restarting calculation. Importing TS from CONTCAR and rewritting MODECAR with contents of NEWMODECAR.")
    TS = read('CONTCAR')
    i = find_last_modecar()
    shutil.copyfile('MODECAR', 'MODECAR%d' % (i+1))
    shutil.copyfile('NEWMODECAR', 'MODECAR')

print("Constraints:")
print("\tc1: Fixing bottom two layers")
del TS.constraints
c1 = FixAtoms(mask=[a.z < 8 for a in TS])
TS.set_constraint(c1)
set_dimer_calc_ZnOCu(atoms_obj = TS)
if testRun:
    run_testRun(TS)
    if not restart:
        run_testRun(TS_vector)
    else:
        TS_vector = TS.copy()
        with open('MODECAR', 'r') as modecar_file:
            for i, line in enumerate(modecar_file):
                TS_vector[i].position += np.array(line)
        run_testRun(TS_vector)

    print("Displaying MODECAR:")
    with open('MODECAR', 'r') as modecar_file:
        print((modecar_file.read()))
else:
    print("Starting dimer calculation.")
    print((TS.get_potential_energy()))
