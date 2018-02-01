import numpy as np
from ase import Atom
from ase.io import read
from ase.io.trajectory import Trajectory
from ase.optimize import LBFGS
from ase.constraints import FixAtoms
from os.path import basename
from run_testRun import run_testRun
from set_calc import set_calc_In2O3
from ase.dft.dos import DOS

testRun = False
file_name_py = basename(__file__)
file_name = file_name_py.replace('.py','')
file_traj = file_name + '.traj'
start_file = '../In2O3_110_clean.traj'

try:
    sys = Trajectory(file_traj, 'r')[-1]
except (IOError, RuntimeError):
    print(("Importing trajectory file from: %s" % start_file))
    sys = read(start_file)
else:
    print(("Importing trajectory file from: %s" % file_traj))

calc = set_calc_In2O3(sys)
print("Constraints:")
print("\tc: Fix bottom two layers.")
c = FixAtoms(mask=[atom.z < 8 for atom in sys])
sys.set_constraint(c)
if testRun == True:
    run_testRun(sys)
else:
    dos = DOS(calc = calc)
    ds = dos.get_dos()
    es = dos.get_energies()
    print("Energy\tDOS")
    for (e, d) in zip(es, ds):
        print(("%f\t%f" % (e, d)))
print(("Completed %s" % file_name_py))    
