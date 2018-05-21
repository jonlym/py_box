# -*- coding: utf-8 -*-
"""
Created on Thu Dec 01 10:51:14 2016

@author: Jonathan Lym
"""
import numpy as np
from scipy.optimize import curve_fit
from py_box3 import interpolate
import py_box3.constants as c
import json

class Shomate(object):
	"""
	Contains the Shomate polynomials and corresponding temperature ranges for a species.

	Attributes:
	-----------
		symbol - string
			Identify the species the object is describing
		T_low - float
			Low temperature (K)
		T_high: float
			High temperature (K)
		phases - dict
			Dictionary of phase objects. Provides higher flexibility for thermodynamics, 
			such as piecewise Shomate polynomials at different temperatures or polymorphs.
			The keys should represent characteristics about the phase (e.g. gamma).
		elements - dict
			Stoichiometric make up of species. Keys should be elements.
			e.g. H2O = {'H': 2, 'O': 1}
	"""
	def __init__(self, symbol, T_low = 0, T_high = 0, a = None, phases = None, elements = None):
		self.symbol = symbol
		if a is None:
			a = np.array(8*[0.])
		if type(phases) is str:
			phase_name = phase
		else:
			phase_name = symbol
		if phases is None:
			phases = {symbol: Phase(symbol = phase_name, T_low = T_low, T_high = T_high, a = a)}
		self.phases = phases	  
		self.elements = elements


	def _get_single_CpoR(self, T, verbose = True, phase = None):
		"""
		Internal function to calculate the dimensionless heat capacity for a single temperature.
		
		Parameters
		----------
			T - float
				Temperature
			verbose - bool
				Whether or not more information should be provided
			phase - string
				Name of phase to calculate the value.
				Other phase options:
					'temperature': Use the shomate polynomial that has the lowest Gibbs energy at T within T_low and T_high [Default]
					'stable': Use the shomate polynomial that has the lowest Gibbs energy at T regardless of T_low and T_high
		Returns
		-------
			CpoR - float
				Dimensionless heat capacity at a single temperature (i.e. Cp/R where R is the molar gas constant)
		"""
		if phase is None or phase == 'temperature':
			phase = self.get_T_phase(T = T)
		elif phase is 'stable':
			phase = self.get_stable_phase(T = T)
		return self.phases[phase].get_CpoR(T = T, verbose = verbose)
	
	def get_CpoR(self, T, verbose = True, phase = None):
		"""
		Calculates the dimensionless heat capacity for a single temperature or a range of temperatures.
		
		Parameters
		----------
			T - float or (N,) ndarray
				Temperature(s)
			verbose - bool
				Whether or not more information should be provided
			phase - string
				Name of phase to calculate the value.
				Other phase options:
					'temperature': Use the shomate polynomial that has the lowest Gibbs energy at T within T_low and T_high [Default]
					'stable': Use the shomate polynomial that has the lowest Gibbs energy at T regardless of T_low and T_high
		Returns
		-------
			CpoR - float or (N,) ndarray
				Dimensionless Heat capacity(ies) corresponding to T (i.e. Cp/R where R is the molar gas constant)
		"""
		try: 
			T_val = iter(T)
		except TypeError:
			#Single value T
			CpoR = self._get_single_CpoR(T = T, verbose = verbose, phase = phase)
		else:
			#List value T
			CpoR = np.array([0.]*len(T))
			for i, T_val in enumerate(T):
				CpoR[i] = self._get_single_CpoR(T = T_val, verbose = verbose, phase = phase)
		return CpoR

	def _get_single_HoRT(self, T, H_correction = False, verbose = True, phase = None):
		"""
		Internal function to calculate the dimensionless enthalpy for a single temperature.
		
		Parameters
		----------
			T - float
				Temperature
			H_correction - bool
				If set to True, subtracts the enthalpy at T = 298 K. Equivalent to H parameter on NIST
			verbose - bool
				Whether or not more information should be provided
			phase - string
				Name of phase to calculate the value.
				Other phase options:
					'temperature': Use the shomate polynomial that has the lowest Gibbs energy at T within T_low and T_high [Default]
					'stable': Use the shomate polynomial that has the lowest Gibbs energy at T regardless of T_low and T_high
		Returns
		-------
			HoRT - float
				Dimensionless enthalpy at a single temperature (i.e. H/RT where R is the molar gas constant)
		"""
		if phase is None or phase == 'temperature':
			phase = self.get_T_phase(T = T)
		elif phase is 'stable':
			phase = self.get_stable_phase(T = T)
		return self.phases[phase].get_HoRT(T = T, verbose = verbose, H_correction = H_correction)

	def get_HoRT(self, T, H_correction = False, verbose = True, phase = None):
		"""
		Calculates the dimensionless enthalpy for a single temperature or a range of temperatures.
		
		Parameters
		----------
			T - float or (N,) ndarray
				Temperature(s)
			H_correction - bool
				If set to True, subtracts the enthalpy at T = 298 K. Equivalent to H parameter on NIST
			verbose - bool
				Whether or not more information should be provided
			phase - string
				Name of phase to calculate the value.
				Other phase options:
					'temperature': Use the shomate polynomial that has the lowest Gibbs energy at T within T_low and T_high [Default]
					'stable': Use the shomate polynomial that has the lowest Gibbs energy at T regardless of T_low and T_high
		Returns
		-------
			HoRT - float or (N,) ndarray
				Dimensionless enthalpy(ies) corresponding to T (i.e. H/RT where R is the molar gas constant)
		"""
		try: 
			T_val = iter(T)
		except:
			#Single value T
			HoRT = self._get_single_HoRT(T = T, H_correction = H_correction, verbose = verbose, phase = phase)
		else:
			#List value T
			HoRT = np.array([0.]*len(T))
			for i, T_val in enumerate(T):
				HoRT[i] = self._get_single_HoRT(T = T_val, H_correction = H_correction, verbose = verbose, phase = phase)
		return HoRT

	def _get_single_SoR(self, T, verbose = True, phase = None):
		"""
		Calculates the dimensionless entropy for a single temperature or a range of temperatures.
		
		Parameters
		----------
			T - float or (N,) ndarray
				Temperature(s)
			verbose - bool
				Whether or not more information should be provided
			phase - string
				Name of phase to calculate the value.
				Other phase options:
					'temperature': Use the shomate polynomial that has the lowest Gibbs energy at T within T_low and T_high [Default]
					'stable': Use the shomate polynomial that has the lowest Gibbs energy at T regardless of T_low and T_high
		Returns
		-------
			SoR - float or (N,) ndarray
				Dimensionless entropy corresponding to T (i.e. S/R where R is the molar gas constant)
		"""		
		if phase is None or phase == 'temperature':
			phase = self.get_T_phase(T = T)
		elif phase is 'stable':
			phase = self.get_stable_phase(T = T)
		return self.phases[phase].get_SoR(T = T, verbose = verbose)

	def get_SoR(self, T, verbose = True, phase = None):
		"""
		Calculates the dimensionless entropy(ies) for a single temperature or a range of temperatures.
		
		Parameters
		----------
			T - float or (N,) ndarray
				Temperature(s)
			verbose - bool
				Whether or not more information should be provided
			phase - string
				Name of phase to calculate the value.
				Other phase options:
					'temperature': Use the shomate polynomial that has the lowest Gibbs energy at T within T_low and T_high [Default]
					'stable': Use the shomate polynomial that has the lowest Gibbs energy at T regardless of T_low and T_high
		Returns
		-------
			SoR - float or (N,) ndarray
				Dimensionless entropy(ies) corresponding to T (i.e. Cp/R where R is the molar gas constant)
		"""		
		try:
			T_val = iter(T)
		except:
			#Single value T
			SoR = self._get_single_SoR(T = T, verbose = verbose, phase = phase)
		else:
			#List value T
			SoR = np.array([0.]*len(T))
			for i, T_val in enumerate(T):
				SoR[i] = self._get_single_SoR(T = T_val, verbose = verbose, phase = phase)
		return SoR

	def get_GoRT(self, T, H_correction = False, verbose = True, phase = None):
		"""
		Calculates the dimensionless Gibbs energy(ies) for a single temperature or a range of temperatures.
		
		Parameters
		----------
			T - float or (N,) ndarray
				Temperature(s)
			H_correction - bool
				If set to True, subtracts the enthalpy at T = 298 K. Equivalent to H parameter on NIST				
			verbose - bool
				Whether or not more information should be provided
			phase - string
				Name of phase to calculate the value.
				Other phase options:
					'temperature': Use the shomate polynomial that has the lowest Gibbs energy at T within T_low and T_high [Default]
					'stable': Use the shomate polynomial that has the lowest Gibbs energy at T regardless of T_low and T_high
		Returns
		-------
			GoRT - float or (N,) ndarray
				Dimensionless Gibbs energy(ies) corresponding to T (i.e. H/RT where R is the molar gas constant)
		"""
		HoRT = self.get_HoRT(T = T, H_correction = H_correction, verbose = verbose, phase = phase)
		SoR = self.get_SoR(T = T, verbose = verbose, phase = phase)
		return HoRT-SoR

	def get_stable_phase(self, T):
		"""
		Use the phase that has the lowest Gibbs energy at T ignoring T_low and T_high.

		Parameters
		----------
			T - float
				Temperature
		Returns
		-------
			min_phase - str
				Name of the phase that has the lowest Gibbs energy at T ignoring T_low and T_high
		"""
		min_phase = None
		min_GoRT = float("inf")
		for name, phase in self.phases.items():
			GoRT = phase.get_GoRT(T = T, verbose = False)
			if GoRT < min_GoRT:
				min_phase = name
				min_GoRT = GoRT
		return min_phase

	def get_T_phase(self, T):
		"""
		Use the phase that has the lowest Gibbs energy at T and within T_low and T_high

		Parameters
		----------
			T - float
				Temperature
		Returns
		-------
			min_phase - str
				Name of the phase that has the lowest Gibbs energy at T ignoring T_low and T_high
		"""
		min_phase = None
		min_GoRT = float("inf")
		for name, phase in self.phases.items():
			GoRT = phase.get_GoRT(T = T, verbose = False)			
			if T <= phase.T_high and T >= phase.T_low and GoRT < min_GoRT:
				min_phase = name
				min_GoRT = GoRT
		return min_phase

	def plot_thermo(self, T_low, T_high, Cp_units = None, S_units = None, H_units = None, G_units = None, phase = None):
		"""
		Plots the heat capacity, enthalpy and entropy in the temperature range specified.
		The units for the plots can be specified by using R

		Parameters
		----------
			T_low - float
				Lower temperature bound
			T_high - float
				Higher temperature bound
			Cp_units - string
				Heat capacity units. Supported units can be found in py_box3.constants.R (but the '/K' has to be omitted)
				e.g. Use 'kJ/mol' to display heat capcity in kiloJoules per mol 
			H_units - string
				Enthalpy units. Supported units can be found in py_box3.constants.R
			S_units - string
				Entropy units. Supported units can be found in py_box3.constants.R (but the '/K' has to be omitted)
				e.g. Use 'kJ/mol' to display entropy in kiloJoules per mol
			G_units - string
				Gibbs free energy units. Supported units can be found in py_box3.constants.R
		Returns
		-------
			fig - matplotlib.figure.Figure object
			ax_Cp - meatplotlib.Axes.axes object
				Axes to heat capacity plot
			ax_H - meatplotlib.Axes.axes object
				Axes to enthalpy plot
			ax_S - meatplotlib.Axes.axes object
				Axes to entropy plot
			ax_G - meatplotlib.Axes.axes object
				Axes to Gibbs energy plot
		"""
		import matplotlib.pyplot as plt

		T = np.linspace(T_low, T_high)
		Cp = self.get_CpoR(T, phase = phase)
		H = self.get_HoRT(T, phase = phase)
		S = self.get_SoR(T, phase = phase)
		G = self.get_GoRT(T, phase = phase)

		#Apply units
		if Cp_units is not None:
			Cp = Cp * c.R('{}/K'.format(Cp_units))
		if S_units is not None:
			S = S * c.R('{}/K'.format(S_units))
		if H_units is not None:
			H = H * c.R(H_units) * T
		if G_units is not None:
			G = G * c.R(G_units) * T

		fig = plt.figure()
		ax_Cp = plt.subplot(411)
		plt.plot(T, Cp, 'r-')
		if units is None:
			plt.ylabel('Cp/R')
		else:
			plt.ylabel('Cp ({}})'.format(Cp_units))
		plt.xlabel('T (K)')
		plt.title('Plots for {} using shomate polynomials.' % self.symbol)

		ax_H = plt.subplot(412)
		plt.plot(T, H, 'g-')
		if units is None:
			plt.ylabel('H/RT')
		else:
			plt.ylabel('H ({})'.format(H_units))
		plt.xlabel('T (K)')

		#Entropy graph
		ax_S = plt.subplot(413)
		plt.plot(T, S, 'b-')
		if units is None:
			plt.ylabel('S/R')
		else:
			plt.ylabel('S ({})'.format(S_units))
		plt.xlabel('T (K)')

		ax_G = plt.subplot(414)
		plt.plot(T, G, 'k-')
		if units is None:
			plt.ylabel('G/RR')
		else:
			plt.ylabel('G (%s)' % G_units)
		plt.xlabel('T (K)')
		return (fig, ax_Cp, ax_H, ax_S, ax_G)

	def get_json_string(self):
		"""
		Converts shomate object to json string

		Returns
		-------
			shomate_json - str
				shomate object in json format
		"""
		shomate_dict = self.__dict__
		for key, phase in self.phases.items():
			phase_dict = phase.__dict__
			phase_dict['a'] = phase.a.tolist()
			shomate_dict['phases'][key] = phase_dict
		return json.dumps(shomate_dict)

	@classmethod
	def fit_shomate_species(cls, symbol, T, Cp, H0, S0, T_ref = c.T0('K'), elements = None):
		"""
		Derives the shomate species from fitting the heat capacity (J/mol/K) and temperature (K) data and including the formation of enthalpy (kJ/mol) and entropy (J/mol/K).

		Paramters
		---------
			symbol - str
				Name of shomate polynomial
			T - (N,) ndarray
				Temperatures (K)
			Cp - (N,) ndarray
				Heat capacities corresponding to T (J/mol/K)
			H0 - float
				Enthalpy at reference temperature (kJ/mol)
			S0 - float
				Entropy at reference temperature (J/mol/K)
			T_ref - float
				Reference temperature (K)
			elements - dict
				Stoichiometric make up of species. Keys should be elements.
				e.g. H2O = {'H': 2, 'O': 1}
		Returns
		-------
			shomate - Shomate object
		"""
		T_low = min(T)
		T_high = max(T)
		t = np.array(T)/1000.
		[a, pcov] = curve_fit(_shomate_Cp, t, np.array(Cp))
		a = np.append(a, [0., 0., 0.])
		a[5] = H0 - _get_HoRT(T = T_ref, a = a)*c.R('kJ/mol/K')*T_ref
		a[6] = S0 - _get_SoR(T = T_ref, a = a)*c.R('J/mol/K')
		a[7] = - _get_HoRT(T = c.T0('K'), a = a)*c.R('kJ/mol/K')*c.T0('K')
		return cls(symbol = symbol, T_low = T_low, T_high = T_high, a = np.array(a), elements = elements)

	@staticmethod
	def read_fund_csv(symbol, csv_path, print_graph = False):
		"""
		Reads a csv file for data that will be fed into function fit_shomate_species.

		Parameters
		----------
			symbol - str
				Name of the shomate polynomial
			csv_path - str
				Path to the csv file containing the values
			print_graph - bool
				Whether or not to show the thermodynamic quantities at the range of temperatures
		Returns
			shomate - Shomate object
		"""
		H0_S0_read = False
		T = []
		Cp = []
		
		print(("Reading from file: %s" % csv_path))
		with open(csv_path, 'r') as csv_file:
			for line in csv_file:
				if line[0] != '!':
					#Split data and remove unnecessary characters
					data = line.split(',')
					T.append(float(data[0]))
					Cp.append(float(data[1]))
					if len(data) > 2 and not H0_S0_read:
						H0 = float(data[2])
						S0 = float(data[3])
						H0_S0_read = True				

		shomate = Shomate.fit_shomate_species(symbol, T, Cp, H0, S0)
		if print_graph:
			T_range = np.linspace(shomate.T_low, shomate.T_high)
			Cp_fit = shomate.get_CpoR(T_range)*c.R('J/mol/K')
			H_fit = shomate.get_HoRT(T_range)*T_range*c.R('kJ/mol/K')
			S_fit = shomate.get_SoR(T_range)*c.R('J/mol/K')
			plt.figure()
			plt.subplot(311)
			plt.plot(T, Cp, 'ro', T_range, Cp_fit, 'b-')
			plt.legend(['NIST Data', 'Fit'])
			plt.xlabel('Temperature (K)')
			plt.ylabel('Cp (J/mol/K)')
			
			plt.subplot(312)
			plt.plot(T_range, H_fit)
			plt.xlabel('Temperature (K)')
			plt.ylabel('H (kJ/mol/K)')

			plt.subplot(313)
			plt.plot(T_range, S_fit)
			plt.xlabel('Temperature (K)')
			plt.ylabel('S (J/mol/K)')
		return shomate
   
def _shomate_Cp(t, A, B, C, D, E):
	"""
	Shomate heat capacity.

	Paramters
	---------
		t - float
			Adjusted temperature (t = T/1000) (K)
		A - float
			Shomate parameter
		B - float
			Shomate parameter
		C - float
			Shomate parameter
		D - float
			Shomate parameter
		E - float
			Shomate parameter
	Returns
	-------
		Cp - float
			Heat capacity at t
	"""
	return A + B*t + C * t**2 + D * t ** 3 + E / t**2 

def generate_Cp_data(Ts, Cps, n = 100):
	"""
	Linearly interpolates between Ts and Cps to create more data. Useful for fitting Cp with few data points

	Parameters
	----------
		Ts - (M,) ndarray
			Input data temperatures
		Cps - (M,) ndarray
			Input heat capacities corresponding to Ts
		n - int
			Number of data points
	Returns
	-------
		Ts_out - (n,) ndarray
			Output temperatures
		Cps_out - (n,) ndarray
			Output heat capacities corresponding to Ts_out
	"""

	Ts_out = np.linspace(min(Ts), max(Ts), n)
	Cps_out = np.zeros(shape = n)

	for i, T in enumerate(Ts_out):
		#Finds nearby values to use for interpolation
		low_i = np.argwhere(T >= Ts)[-1]
		high_i = np.argwhere(T <= Ts)[0]
		low_Cp = Cps[low_i]
		high_Cp = Cps[high_i]
		low_T = Ts[low_i]
		high_T = Ts[high_i]

		Cps_out[i] = interpolate(x_low = low_T, x_high = high_T, y_low = low_Cp, y_high = high_Cp, x = T)
	return (Ts_out, Cps_out)

def _get_HoRT(T, a, H_correction = False):
	"""
	Internal function to calculate dimensionless enthalpy.
	
	Parameters
	----------
		T - float or (N,) ndarray
			Temperature(s)
		a - (8,) ndarray
			Shomate parameters
		H_correction - bool
			If set to True, subtracts the enthalpy at T = 298 K. Equivalent to H parameter on NIST
	Returns
	-------
		HoRT - float or (N,) ndarray
			Dimensionless enthalpy corresponding to T (i.e. H/RT where R is the molar gas constant)
	"""
	t = T/1000.
	T_arr = np.array([t, t ** 2 / 2, t ** 3 / 3, t ** 4 / 4, -1/t, 1., 0., 1.])
	return np.dot(T_arr, a)/(c.R('kJ/mol/K')*T)

def _get_SoR(T, a):
	"""
	Internal function to calculate dimensionless entropy.
	
	Parameters
	----------
		T - float or (N,) ndarray
			Temperature(s)
		a - (8,) ndarray
			Shomate parameters
	Returns
	-------
		HoRT - float or (N,) ndarray
			Dimensionless entropy corresponding to T (i.e. S/R where R is the molar gas constant)
	"""
	t = T/1000.		
	T_arr = np.array([np.log(t), t, t ** 2 / 2., t ** 3 / 3., -0.5 * ( 1 / t ) ** 2, 0., 1., 0.])
	return np.dot(T_arr, a)/c.R('J/mol/K')

class Phase(object):
	"""
	Contains the Shomate polynomial for a phase. This may include a phase transition,
	a polymorph, an empirical piece-wise polynomial, etc.

	Attributes
	----------
		symbol - str
			Name of phase used to identify it
		T_low - float
			Low Temperature (K)
		T_high - float
			High Temperature (K)
		a - (8,) ndarray
			Shomate polynomial
	"""

	def __init__(self, symbol = '', T_low = 0., T_high = 0., a = None):
		self.symbol = symbol
		self.T_low = T_low
		self.T_high = T_high
		self.a = a

	def get_CpoR(self, T, verbose = True):
		"""
		Calculate the dimensionless heat capacity for a single temperature.
		
		Parameters
		----------
			T - float
				Temperature
			verbose - bool
				Whether or not more information should be provided
		Returns
		-------
			CpoR - float
				Dimensionless heat capacity at a single temperature (i.e. Cp/R where R is the molar gas constant)
		"""
		t = T/1000.		
		T_arr = np.array([1, t, t ** 2, t ** 3, (1/t) ** 2, 0, 0, 0])
		if verbose:
			if T < self.T_low:
				print(("Warning. Input temperature (%f) lower than T_low (%f) for species %s" % (T, self.T_low, self.symbol)))
			elif T > self.T_high:
				print(("Warning. Input temperature (%f) higher than T_high (%f) for species %s" % (T, self.T_high, self.symbol)))
		return np.dot(T_arr, self.a)/c.R('J/mol/K')
	
	def get_HoRT(self, T, H_correction = False, verbose = True):
		"""
		Calculate the dimensionless enthalpy for a single temperature.
		
		Parameters
		----------
			T - float
				Temperature
			H_correction - bool
				If set to True, subtracts the enthalpy at T = 298 K. Equivalent to H parameter on NIST
			verbose - bool
				Whether or not more information should be provided
		Returns
		-------
			HoRT - float
				Dimensionless enthalpy at a single temperature (i.e. H/RT where R is the molar gas constant)
		"""		
		t = T/1000.
		if H_correction:
			T_arr = np.array([t, t ** 2 / 2, t ** 3 / 3, t ** 4 / 4, -1/t, 1., 0., 1.])
		else:			
			T_arr = np.array([t, t ** 2 / 2, t ** 3 / 3, t ** 4 / 4, -1/t, 1., 0., 0.])
		if verbose:
			if T < self.T_low:
				print(("Warning. Input temperature (%f) lower than T_low (%f) for species %s" % (T, self.T_low, self.symbol)))
			elif T > self.T_high:
				print(("Warning. Input temperature (%f) higher than T_high (%f) for species %s" % (T, self.T_high, self.symbol)))
		return np.dot(T_arr, self.a)/(c.R('kJ/mol/K')*T)

	def get_SoR(self, T, verbose = True):
		"""
		Calculates the dimensionless entropy for a single temperature.
		
		Parameters
		----------
			T - float or (N,) ndarray
				Temperature(s)
			verbose - bool
				Whether or not more information should be provided
		Returns
		-------
			SoR - float
				Dimensionless entropy corresponding to T (i.e. S/R where R is the molar gas constant)
		""" 
		t = T/1000.		
		T_arr = np.array([np.log(t), t, t ** 2 / 2., t ** 3 / 3., -0.5 * ( 1 / t ) ** 2, 0., 1., 0.])
		if verbose:
			if T < self.T_low:
				print(("Warning. Input temperature (%f) lower than T_low (%f) for species %s" % (T, self.T_low, self.symbol)))
			elif T > self.T_high:
				print(("Warning. Input temperature (%f) higher than T_high (%f) for species %s" % (T, self.T_high, self.symbol)))
		return np.dot(T_arr, self.a)/c.R('J/mol/K')

	def get_GoRT(self, T, H_correction = False, verbose = True):
		"""
		Calculate the dimensionless Gibbs free energy for a single temperature.
		
		Parameters
		----------
			T - float
				Temperature
			H_correction - bool
				If set to True, subtracts the enthalpy at T = 298 K. Equivalent to H parameter on NIST
			verbose - bool
				Whether or not more information should be provided
		Returns
		-------
			GoRT - float
				Dimensionless enthalpy at a single temperature (i.e. H/RT where R is the molar gas constant)
		"""
		HoRT = self.get_HoRT(T, H_correction = H_correction, verbose = verbose)
		SoR = self.get_SoR(T, verbose)
		return HoRT-SoR