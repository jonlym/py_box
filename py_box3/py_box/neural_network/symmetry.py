import numpy as np

def get_G2(i, atoms, eta, gamma, xi, R_c, R_s):
	G2 = 0.
	for atom in atoms:
		if atom.index == i:
			continue

		distance = atoms.get_distance(atom.index, i, mic = True)
		f = get_f(distance, R)
		G2 += np.exp(-1. * eta * (distance - R_s)**2) * f
	return G2

def get_G4(i, atoms, eta, gamma, xi, R):
	G4 = 0.
	for atom_j in atoms:
		for atom_k in atoms:
			#Check if any of the indices are the same
			if len(set((i, atom_j.index, atom_k.index))) != 3:
				continue

			distance_ij = atoms.get_distance(atom_j.index, i, mic = True)
			distance_ik = atoms.get_distance(atom_k.index, i, mic = True)
			f_ij = get_f(distance_ij, R)
			f_ik = get_f(distance_ik, R)
			cos_theta = get_cos_theta(atoms, i, atom_j.index, atom_k.index)

			G4 += (1. + gamma * cos_theta)**xi * np.exp(-1. * eta * (distance_ij**2 + distance_ik**2)) * f_ij * f_ik
	return 2.**(1.-xi) * G4

def get_f(distance, R):
	if distance <= R:
		f = 0.5 * (np.cos(np.pi * distance / R) + 1.)
	else:
		f = 0.
	return f

def get_cos_theta(atoms, i, j, k):
	r_ij = atoms.get_distance(i, j, mic = True, vector = True)
	r_ik = atoms.get_distance(i, k, mic = True, vector = True)
	return np.vdot(r_ij, r_ik) / (np.linalg.norm(r_ij) * np.linalg.norm(r_ik))