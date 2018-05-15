import collections

class Reactions(collections.UserList):
    """
    A User Dictionary object that stores Reaction objects. The key is the symbol name.
    """	

    def get_stoichiometric_matrix(thermos):
    	"""
		Returns the stoichiometric matrix.
		Parameters
		----------
		thermos - Mapping object of thermodynamic objects (e.g. Shomates)
			Container that holds the thermodynamic objects. The container should use the
			species found in stoichiometry as keys

		Returns
		-------
			stoich_mat - (n_reactions, n_species) ndarray
			Stoichiometric matrix
		"""
    	species = list(thermos.keys())
    	stoich_mat = np.zeros((len(reactions), len(thermos)))
    	for i, reaction in enumerate(self):
    		for symbol, coefficient in reaction.items():
    			j = species.index(symbol)
    			stoich_mat[i, j] = coefficient
    	return stoich_mat