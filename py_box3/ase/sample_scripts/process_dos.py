import numpy as np

from ase.io import read
from pykit.densityofstates import dosread, dosplot

# dosread returns a DensityOfStates object.  Pass the DOSCAR and either the
# POSCAR or CONTCAR file
dos = dosread('DOSCAR', posfile='CONTCAR')
print((dos.dos_data))
# Grab the occupied p energy for every single Al
print(('E_In(d states) = %.2f\n' % dos.get_band_energy('d', atomtype='In', erange = [-30, 0])))

# Now, here is where things get a bit weird, 'cause I apparently cannot code 
# up an easy to use interface :) 

# Lets grab the p energy of a specific surface atom.  Here, I'll randomly pick
# Al27.  Now, for whatever reason that has been lost to the depths of time, I 
# opted to reference each atom by an atom label and index (i.e. 'Al27').

# The DensityOfStates object can retrieve the energy via the get_band_energy
# class function.  By default, it will grab all the atoms of a particular type.
# This behavior can be modified via the exclude/include options, which either
# can be a string or a list.

# So to grab a single atom, we will need to construct an exclude list in this
# fashion:

atoms = read('CONTCAR')
exclude = list()
all_atoms = list()
for atom in atoms:
    all_atoms.append('%s%i' % (atom.symbol, atom.index))
    if atom.index != 24:
        exclude.append('%s%i' % (atom.symbol, atom.index))
print(all_atoms)
print(('E_In(d states) = %.2f\n' % dos.get_band_energy('d', atomtype='In', include = all_atoms)))

# Do the integration!  Since we want both the valence and conduction regions, we set
# some integration limits (this will be case dependent!
print(('E_s = %.2f' % dos.get_band_energy('s', atomtype='In', include='In24', exclude=exclude, erange=[-30, 0.])))
print(('E_p = %.2f' % dos.get_band_energy('p', atomtype='In', include='In24', exclude=exclude, erange=[-30, 0.])))
print(('E_d = %.2f' % dos.get_band_energy('d', atomtype='In', include='In24', exclude=exclude, erange=[-30, 0.])))
print(('E^*_s = %.2f' % dos.get_band_energy('s', atomtype='In', include='In24', exclude=exclude, erange=[0, 20.])))
print(('E^*_p = %.2f' % dos.get_band_energy('p', atomtype='In', include='In24', exclude=exclude, erange=[0, 20.])))
print(('E^*_d = %.2f' % dos.get_band_energy('d', atomtype='In', include='In24', exclude=exclude, erange=[0, 20.])))

# Compare these values to Table S3 in the SI from: J. Phys. Chem. C, 118, 12899 (2014)


###############################################################################

# Now, lets plot up some DOS
# First get the energy grid.  This should be zero'ed to the Fermi level
energies = dos.get_energy_grid()

# get_atom_dos will return the DOS for a specific atom label as a numpy array
# Each column corresponds to the s, p, and d orbitals
d_in24 = dos.get_atom_dos('In24')

# Pretty self-explanatory.  This function calls matplotlib, so if you don't
# have or don't want to use, ignore
#dosplot('site3.eps', energies=energies, dos=[d_al27[:, 0], d_al27[:, 1], d_al27[:, 2]], legend=['s', 'p', 'd'])

# If you want to extract out the DOS, and then read them into something else, 
# this shows you how to save both the energies and DOS as a text file.
shape = d_in24.shape

data = np.empty((shape[0], shape[1] + 1))
data[:, 0] = energies
data[:, 1:] = d_in24
np.savetxt('dos_output.txt', data, fmt='%15.5f')



