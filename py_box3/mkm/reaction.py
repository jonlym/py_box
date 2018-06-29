import numpy as np
from py_box3 import any_alpha
from py_box3 import constants as c

class Reaction(object):
	"""
	Contains a chemical reaction.
	Attributes
	----------
		stoichiometry: dict
			Stoichiometry of the reaction. Negative values indicate reactants and positive values indicate products. 
			The keys should be the elements located in a thermodynamic dictionary, such as a Shomates object
		A: float
			Pre-exponential factor in s. Default value is kb/h
		EoRT: float
			Dimensionless activation energy (T = 298 K)
		BEP: (2,) ndarray
			Bronsted Evans Polyani polynomial coefficients. EoRT = BEP[0]*delta_HoRT + BEP[1]
		beta: float
			Power to be used for temperature correction. i.e. (T/T_ref)^beta
		T_ref: float
			Reference temperature to be used for temperature correction. i.e. (T/T_ref)^beta
		use_sticking_coefficient: boolean
			Calculate A using the sticking coefficient. Useful for adsorption reactions.
		s: float
			Sticking coefficient for adsorption reactions		
	"""
	def __init__(self, stoichiometry, A = c.kb('J/K')/c.h('J s'), EoRT = None, BEP = None, beta = 0., T_ref = 1., use_sticking_coefficient = False, s = 0.5):
		self.stoichiometry = stoichiometry
		self.A = A
		self.EoRT = EoRT
		self.BEP = BEP
		self.beta = beta
		self.T_ref = T_ref
		self.use_sticking_coefficient = use_sticking_coefficient
		self.s = s

	def get_A_from_sticking_coefficient(T, thermos):

		#Find gas-phase reactant
		for species, coefficient in self.stoichiometry.items():
			if coefficient < 0 and thermos[species].is_gas:
				#Calculate molecular weight
				MW = c.get_molecular_weight(thermos[species].elements)
				A = self.s * np.sqrt(c.R('J/mol/K')*T/(2*np.pi*MW))
				return A


	def get_delta_HoRT(self, T, thermos):
		"""
		Calculates the enthalpy of reaction at a certain temperature.
		Arguments
		---------
			T - float or (N,) ndarray
				Temperature (K)
			thermos - dict-like object of thermodynamic objects (e.g. Shomates)
				Container that holds the thermodynamic objects. The container should use the
				species found in stoichiometry as keys and the objects should have the method: get_HoRT()
		
		Returns
		-------
			HoRT - float or (N,) ndarray
				Dimensionless enthalpy of reaction
		"""

		#Initialize HoRT
		try:
			iter(T)
		except TypeError:
			HoRT = 0.
		else:
			HoRT = np.zeros(len(T))

		for species, coefficient in self.stoichiometry.items():
			HoRT += coefficient * thermos[species].get_HoRT(T = T)
		return HoRT

	def get_delta_SoR(self, T, thermos):
		"""
		Calculates the entropy of reaction at a certain temperature.
		Arguments
		---------
			T - float or (N,) ndarray
				Temperature (K)
			thermos - dict-like object of thermodynamic objects (e.g. Shomates)
				The container should use the keys of stoichiometry as keys and the objects should have the method: get_SoR()
		
		Returns
		-------
			HoRT - float or (N,) ndarray
				Dimensionless enthalpy of reaction
		"""

		#Initialize HoRT
		try:
			iter(T)
		except TypeError:
			SoR = 0.
		else:
			SoR = np.zeros(len(T))

		for species, coefficient in self.stoichiometry.items():
			SoR += coefficient * thermos[species].get_SoR(T = T)
		return SoR


	def get_delta_GoRT(self, T, thermos):
		"""
		Calculates the Gibbs energy of reaction at a certain temperature.
		Arguments
		---------
			T - float or (N,) ndarray
				Temperature (K)
			thermos - dict-like object of thermodynamic objects (e.g. Shomates)
				Container that holds the thermodynamic objects. The container should use the
				species found in stoichiometry as keys and the objects should have the method: get_GoRT()
		
		Returns
		-------
			GoRT - float or (N,) ndarray
				Dimensionless Gibbs energy of reaction
		"""

		#Initialize HoRT
		try:
			iter(T)
		except TypeError:
			GoRT = 0.
		else:
			GoRT = np.zeros(len(T))

		for species, coefficient in self.stoichiometry.items():
			GoRT += coefficient * thermos[species].get_GoRT(T = T)
		return GoRT

	def get_EoRT_from_BEP(self, T, thermos):
		"""
		Calculates the activation energy from BEP relationship.
		Arguments
		---------
			T - float or (N,) ndarray
				Temperature (K)
			thermos - dict-like object of thermodynamic objects (e.g. Shomates)
				Container that holds the thermodynamic objects. The container should use the
				species found in stoichiometry as keys and the objects should have the method: get_HoRT()
		Returns
		-------
			EoRT - float or (N,) ndarray
				Dimensionless activation energy
		"""
		HoRT = self.get_HoRT(T = T, thermos = thermos)
		return np.polyval(self.BEP, HoRT)

	def get_forward_rate_constant(self, T, thermos = None, use_BEP = False):
		"""
		Calculates the forward rate constant, k_fwd
		Arguments
		---------
			T - float
				Temperature (K)
			thermos - dict-like object of thermodynamic objects (e.g. Shomates)
				Container that holds the thermodynamic objects. The container should use the
				species found in stoichiometry as keys and the objects should have the method: get_HoRT()
			use_BEP - boolean
				Whether to use BEP relationship
		Returns
		-------
			k_fwd - float
				Forward rate constant
		"""
		if use_BEP:
			EoRT = self.get_EoRT_from_BEP(T = T, thermos = thermos)
		else:
			EoRT = self.EoRT

		if use_sticking_coefficient:
			A_fwd = self.get_A_from_sticking_coefficient(T = T)
		else:
			A_fwd = self.A

		return self.A * (T/self.T_ref)**self.beta * np.exp(EoRT * c.T0('K') / T)

	def get_reverse_rate_constant(self, T, thermos, use_BEP = False):
		"""
		Calculates the reverse rate constant, k_rev
		Arguments
		---------
			T - float
				Temperature (K)
			thermos - dict-like object of thermodynamic objects (e.g. Shomates)
				Container that holds the thermodynamic objects. The container should use the
				species found in stoichiometry as keys and the objects should have the methods: get_HoRT() and get_SoRT()
			use_BEP - boolean
				Whether to use BEP relationship
		Returns
		-------
			k_rev - float
				Reverse rate constant
		"""

		#Get activation energy of reverse direction
		if use_BEP:
			EoRT = self.get_EoRT_from_BEP(T = T, thermos = thermos) - self.get_HoRT(T = T, thermos = thermos)
		else:
			EoRT = self.EoRT - self.get_HoRT(T = T, thermos = thermos)

		#Get preexponential factor of reverse direction
		if use_sticking_coefficient:
			A_fwd = self.get_A_from_sticking_coefficient(T = T)
		else:
			A_fwd = self.A
		A_rev = A_fwd*np.exp(self.get_SoR(T = T, thermos = thermos))

		return A_rev * (T/self.T_ref)**self.beta * np.exp(EoRT * c.T0('K') / T)

	def get_forward_rate(self, T, concentrations, thermos = None, use_BEP = False):
		"""
		Calculates the forward rate, r_fwd
		Arguments
		---------
			T - float
				Temperature (K)
			concentrations - dict-like object
				Holds the concentrations of each species. The keys should be the same as used in thermos
			thermos - dict-like object of thermodynamic objects (e.g. Shomates)
				Container that holds the thermodynamic objects. The container should use the
				species found in stoichiometry as keys and the objects should have the method: get_HoRT()
			use_BEP - boolean
				Whether to use BEP relationship
		Returns
		-------
			r_fwd - float
				Forward reaction rate
		"""
		k_fwd = self.get_forward_rate_constant(T = T, thermos = thermos, use_BEP = use_BEP)
		r_fwd = k_fwd
		for symbol, coefficient in self:
			if coefficient > 0:
				r_fwd *= concentrations[symbol] ** coefficient
		return r_fwd

	def get_reverse_rate(self, T, concentrations, thermos = None, use_BEP = False):
		"""
		Calculates the reverse rate, r_rev
		Arguments
		---------
			T - float
				Temperature (K)
			concentrations - dict-like object
				Holds the concentrations of each species. The keys should be the same as used in thermos
			thermos - dict-like object of thermodynamic objects (e.g. Shomates)
				Container that holds the thermodynamic objects. The container should use the
				species found in stoichiometry as keys and the objects should have the methods: get_HoRT() and get_SoR()
			use_BEP - boolean
				Whether to use BEP relationship
		Returns
		-------
			r_rev - float
				Reverse reaction rate
		"""
		k_rev = self.get_reverse_rate_constant(T = T, thermos = thermos, use_BEP = use_BEP)
		r_rev = k_rev
		for symbol, coefficient in self:
			if coefficient < 0:
				r_rev *= concentrations[symbol] ** np.abs(coefficient)
		return r_rev


	def get_net_rate(self, T, concentrations, thermos, use_BEP = False):
		"""
		Calculates the reverse rate, r_rev
		Arguments
		---------
			T - float
				Temperature (K)
			concentrations - dict-like object
				Holds the concentrations of each species. The keys should be the same as used in thermos
			thermos - dict-like object of thermodynamic objects (e.g. Shomates)
				Container that holds the thermodynamic objects. The container should use the
				species found in stoichiometry as keys and the objects should have the methods: get_HoRT() and get_SoR()
			use_BEP - boolean
				Whether to use BEP relationship
		Returns
		-------
			r_net - float
				Net reaction rate
		"""
		r_fwd = self.get_forward_rate(T = T, concentrations = concentrations, thermos = thermos, use_BEP = use_BEP)
		r_rev = self.get_reverse_rate(T = T, concentrations = concentrations, thermos = thermos, use_BEP = use_BEP)
		return r_fwd - r_rev

	def get_pe_ratio(self, T, concentrations, thermos, use_BEP = False):
		"""
		Calculates the partial equilibrium index
		Arguments
		---------
			T - float
				Temperature (K)
			concentrations - dict-like object
				Holds the concentrations of each species. The keys should be the same as used in thermos
			thermos - dict-like object of thermodynamic objects (e.g. Shomates)
				Container that holds the thermodynamic objects. The container should use the
				species found in stoichiometry as keys and the objects should have the methods: get_HoRT() and get_SoR()
			use_BEP - boolean
				Whether to use BEP relationship
		Returns
		-------
			pe_ratio - float
				Partial equilibrium index 
		"""
		r_fwd = self.get_forward_rate(T = T, concentrations = concentrations, thermos = thermos, use_BEP = use_BEP)
		r_rev = self.get_reverse_rate(T = T, concentrations = concentrations, thermos = thermos, use_BEP = use_BEP)
		return r_fwd/(r_fwd + r_rev)

	def check_stoichiometry(self, thermos):
		"""
		Checks the stoichiometry of a reaction. Will throw a ValueError if the reaction is not balanced.

		Arguments
		---------
			thermos - dict-like object of thermodynamic objects (e.g. Shomates)
				Container that holds the thermodynamic objects. The container should use the
				species found in stoichiometry as keys and the objects should have the methods: get_HoRT() and get_SoR()
		Raises
		-------
			ValueError
				If the reaction is not balanced

		"""
		balance = {}
		for species, coefficient_species in self.coefficients.items():
			for element, coefficient_element in thermos[species].elements.items():
				try:
					balance[species] += coefficient_element * coefficient_species
				except KeyError:
					balance[species] = coefficient_element * coefficient_species
		if any(value != 0 for key, value in balance.items()):
			raise ValueError('Equation not balanced.')

def parse_reaction(reaction):
	"""
	Parses a reaction
	Args
		reaction - str
			Balanced chemical reaction
	Returns
		dict
			Reaction translated to dictionary where the keys are the species and the values is the stoichoimetric coefficient
	"""
	reaction_dict = {}
	reaction_side = -1 #-1 for reactants, +1 for products
	separated_reaction = reaction.split(' ')
	for i, field in enumerate(separated_reaction):
		if any_alpha(field):
			if i == 0:
				reaction_dict[field] = 1. * reaction_side
			else:
				try:
					reaction_dict[field] = float(separated_reaction[i-1]) * reaction_side
				except ValueError:
					reaction_dict[field] = 1. * reaction_side
		elif field == '->':
			reaction_side = 1.
	return reaction_dict

def write_reaction(reaction_dict):
	"""
	Converts the reaction elements into its stoichiometric form
	Args
		reaction_dict - dict
			Reaction as a dictionary where the keys are the species and the values are the stoichiometric coefficient
	Returns
		str
			Chemical reaction as a string
	"""
	reaction = ''
	#Add reactants
	for specie, coefficient in reaction_dict.items():
		if coefficient < 0.:
			reaction = '{}{:.3f}{} + '.format(reaction, -coefficient, specie)
	else:
		reaction = reaction[:-2]

	reaction = '{}-> '.format(reaction)

	for specie, coefficient in reaction_dict.items():
		if coefficient > 0.:
			reaction = '{}{:.3f}{} + '.format(reaction, coefficient, specie)
	else:
		reaction = reaction[:-2]
	return reaction
