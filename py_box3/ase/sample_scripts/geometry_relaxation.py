import numpy as np
from ase import Atom
from ase.build import surface
from ase.build import sort
from ase.constraints import FixAtoms
from ase.io import read
from ase.io import write
from ase.io.trajectory import Trajectory
from ase.optimize import LBFGS
from os.path import basename
from run_testRun import run_testRun
from set_calc import set_calc_In2O3

testRun = False
file_name_py = basename(__file__)
file_name = file_name_py.replace('.py','')
file_traj = file_name + '.traj'
start_file = '../In2O3_110_clean.traj'

mode = 'a'
try:
    sys = Trajectory(file_traj, 'r')[-1]
except (IOError, RuntimeError):
    print(("Importing trajectory file from: %s" % start_file))
    sys = read(start_file)
    sys.center(vacuum = 5.0, axis = 2)
else:
    print(("Importing trajectory file from: %s" % file_traj))

set_calc_In2O3(sys)
print("Constraints:")
print("	c1: Fix bottom two layers.")
avg_z = np.mean([atom.z for atom in sys])
c1 = FixAtoms(mask=[atom.z < avg_z for atom in sys])
sys.set_constraint(c1)

if testRun == True:
    run_testRun(sys)
else:
    geo_traj = Trajectory(file_traj, mode, sys)
    dyn = LBFGS(sys, trajectory = geo_traj)
    dyn.run(fmax = 0.05)
    energy = sys.get_potential_energy()
    print(("Energy: %f" % energy))
print(("Completed %s" % file_name_py))
