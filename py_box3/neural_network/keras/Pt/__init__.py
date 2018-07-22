import jenkspy
import numpy as np
from copy import copy
from py_box3 import get_unique_list, base10_to_basen

def get_n_layers(atoms, n_metal_layers = 5):
	"""Returns the number of layers in the system."""
	atom_symbols = atoms.get_chemical_symbols()
	if 'N' in atom_symbols:
		n_layers = n_metal_layers + 2
	elif 'O' in atom_symbols:
		n_layers = n_metal_layers + 1
	else:
		n_layers = n_metal_layers
	return n_layers

def get_breaks(atoms, n_layers):	
	atom_zs = np.array([atom.z for atom in atoms])
	try:
		breaks = jenkspy.jenks_breaks(atom_zs, nb_class = n_layers)
	except ValueError:
		breaks = atom_zs
	return get_unique_list(sorted(breaks))

def assign_layers(atoms, breaks = None, spacing = 2.302, n_metal_layers = 5):
	if breaks is None:
		breaks = simple_get_breaks(atoms = atoms, spacing = spacing, n_metal_layers = n_metal_layers)
	return [np.argmin(np.abs(atom.z - breaks)) for atom in atoms]
	# layers = []
	# for atom in atoms:
	# 	for i, br in enumerate(breaks):			
	# 		if atom.z <= br:
	# 			if i == 0:
	# 				layers.append(i)
	# 			else:
	# 				layers.append(i-1)
	# 			break
			# #If the z coordinate is the upper break, need to put in lower class
			# elif (atom.z == br) and (i == len(breaks)-1):
			# 	layers.append(i - 1)
			# 	break
	# return layers
	#return [np.where(atom.z < breaks)[0][0]-1 for atom in atoms]

	# if len(breaks) != n_layers:
	# 	breaks = breaks[:-1]


def del_lower_layers(atoms, spacing = 2.302, n_metal_layers = 5, n_keep = 1):
	# n_layers = get_n_layers(atoms = atoms, n_metal_layers = n_metal_layers)
	# breaks = get_breaks(atoms = atoms, n_layers = n_layers)
	# layers = assign_layers(atoms = atoms, breaks = breaks)
	layers = assign_layers(atoms = atoms, spacing = spacing)

	del_indices = [atom.index for atom, layer in zip(atoms, layers) if layer < n_metal_layers - n_keep]
	del atoms[del_indices]

def get_fcc_positions(gcn, site_type, spacing = 2.302, n_metal_layers = 2, metals = ['Pt'], offset = 0.):
	atoms = gcn.atoms
	site_occupancies = []
	# n_layers = get_n_layers(atoms = atoms, n_metal_layers = n_metal_layers)
	# breaks = get_breaks(atoms = gcn.atoms, n_layers = n_layers)
	# layers = assign_layers(atoms = atoms, breaks = breaks)

	layers = assign_layers(atoms = atoms, spacing = spacing)

	fcc_positions = []
	for i, (atom_i, layer_i, neighbors_i) in enumerate(zip(atoms, layers, gcn.neighbors)):
		#Atom must be the metal or on the top layer
		if atom_i.symbol not in metals or layer_i != n_metal_layers - 1:
			continue

		neighbors_i = set(neighbors_i)
		for j in neighbors_i:
			#Atom must be the metal or on the top layer
			if atoms[j].symbol not in metals or layers[j] != n_metal_layers - 1:
				continue

			neighbors_j = set(gcn.neighbors[j])
			#Find the common neighbors
			neighbors_ij = set(neighbors_i).intersection(neighbors_j)
			#Move onto the next neighbor if they do not share neighbors
			if len(neighbors_ij) == 0:
				continue

			for k in neighbors_ij:
				#Move to next shared neighbor if not on the surface
				if layers[k] != n_metal_layers - 1:
					continue

				neighbors_k = set(gcn.neighbors[k])
				neighbors_ijk = neighbors_ij.intersection(neighbors_k)
				#None of their neighbors can be in the second row			
				if len(neighbors_ijk) > 0 and any([layers[l] == n_metal_layers - 2 for l in neighbors_ijk]):
					continue

				fcc_position = np.mean([atoms[i].position, atoms[j].position, atoms[k].position], axis = 0) + np.array([0., 0., offset])

				#Skip dupliates
				if any([np.array_equal(fcc_position, x) for x in fcc_positions]):
					continue

				fcc_positions.append(fcc_position)
				#See if the site is occupied
				for m in neighbors_ijk:
					if atoms[m].symbol in site_type.keys():
						site_occupancies.append(site_type[atoms[m].symbol])
						break
				else:
					site_occupancies.append(0)

	return (fcc_positions, site_occupancies)

def simple_get_breaks(atoms, spacing = 2.302, n_metal_layers = 5):
	z_min = np.min([atom.z for atom in atoms])
	breaks = z_min + np.array([range(n_metal_layers)]) * spacing
	return breaks

def simple_get_fcc_sites(atoms, site_type, n_metal_layers, spacing = 2.302, metals = ['Pt'], offset = 0., tol = {'O': 0.7, 'N': 0.7}, n_dimensions = 3):
	ads_positions = {}
	for element in site_type.keys():
		ads_positions[element] = [atom.position for atom in atoms if atom.symbol == element]

	fcc_occupancies = []
	fcc_positions = []
	fcc_periodic_positions = []

	layers = assign_layers(atoms = atoms, spacing = spacing, n_metal_layers = n_metal_layers)
	for atom, layer in zip(atoms, layers):
		#Atom must be on the third layer
		if layer == n_metal_layers - 3:
			#Record position
			fcc_position = atom.position + np.array([0., 0., spacing * 2 + offset])
			fcc_positions.append(fcc_position)

			#Check periodic images for site occupancy
			periodic_positions = get_periodic_positions(atoms = atoms, position = fcc_position, n_dimensions = n_dimensions)
			#fcc_periodic_positions.append(periodic_positions)

			ads_distances = {}
			for periodic_position in periodic_positions:
				occupancy_found = False
				for element, position in ads_positions.items():
					ads_distances[element] = np.array([np.linalg.norm(position - periodic_position, ord = 2) for position in ads_positions[element]])

				for element, distances in ads_distances.items():
					if np.any(distances <= tol[element]):
						fcc_occupancies.append(site_type[element])
						occupancy_found = True
						break

				if occupancy_found:
					break
			else:
				fcc_occupancies.append(-1)

	return (fcc_positions, fcc_occupancies)


def get_periodic_positions(position, atoms = None, cell = None, n_dimensions = 3):
	if atoms is not None:
		cell = atoms.get_cell()
	offsets = np.array([base10_to_basen(num = i, n = 3, width = n_dimensions)-1 for i in range(3**n_dimensions)])

	periodic_positions = []
	for offset in offsets:
		periodic_position = np.dot(cell.T, offset) + position
		periodic_positions.append(periodic_position)

	return np.array(periodic_positions)