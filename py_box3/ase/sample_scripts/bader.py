
import numpy as np
from ase.io import read
from ase.io.trajectory import Trajectory
from ase.optimize import LBFGS
from ase.constraints import FixAtoms
from os.path import basename
from os import system
from run_testRun import run_testRun
from set_calc import set_bader_calc

testRun = False
file_name_py = basename(__file__)
file_name = file_name_py.replace('.py','')
file_traj = file_name + '.traj'
start_file = '../In2O3_110_H-1.traj'
mode = 'a'
try:
    sys = Trajectory(file_traj, 'r')[-1]
except (IOError, RuntimeError):
    print(("Importing trajectory file from: %s" % start_file))
    sys = read(start_file)
    mode = 'w'
else:
    print(("Importing trajectory file from: %s" % file_traj))

set_bader_calc(sys, sigma = 0.05)
print("Constraints:")
print("	c: Fix bottom two layers.")
c = FixAtoms(mask=[atom.z < 8 for atom in sys])
sys.set_constraint(c)
if testRun == True:
    run_testRun(sys)
else:
    geo_traj = Trajectory(file_traj, mode, sys)
    dyn = LBFGS(sys, trajectory = geo_traj)
    dyn.run(fmax = 0.05)
    energy = sys.get_potential_energy()
    print(("Energy: %f" % energy))
    print('Generating CHGCAR_sum for Bader analysis...')
    system('chgsum.pl AECCAR0 AECCAR2')
    print('Performing Bader analysis (bader CHGCAR -ref CHGCAR_sum)')
    system('bader CHGCAR -ref CHGCAR_sum')
print(('Completed %s' % file_name_py))
