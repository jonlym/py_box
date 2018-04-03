# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 16:38:29 2016

@author: Jonathan Lym
"""
from py_box3 import get_unique_list

calc_dict = {
	'ZnOCu': {
		'xc': "PBE",
		'kpts': (3,3,1), 
		'encut': 400,
		'ismear': 0,
		'sigma': 0.1,
		'ediff': 1e-4,
		'prec': 'normal',
		'lcharg': False,
		'lwave': True,
		'nelmin': 4,
		'nelmdl': 6,
		'npar': 2,
		'algo': 'fast',
		'lreal': 'auto',
		'ispin': 2,
		'magmom': 0,
		'istart': 0,
		},
	'In2O3': {
		'xc': "PBE",
		'setups': {'In': '_d'},
		'kpts': (4,3,1),
		'encut': 400,
		'ismear': 0,
		'sigma': 0.05,
		'ediff': 1e-4,
		'prec': 'normal',
		'lcharg': False,
		'lwave': True,
		'nelmin': 4,
		'nelmdl': 6,
		'npar': 2,
		'algo': 'fast',
		'lreal': 'auto',
		'ispin': 2,
		'gamma': True,
		'istart': 0,
		'nelm': 400,
		'nelmdl': -10,
		},
	'ZnOCu_bader': {
		'xc': "PBE",
		'kpts': (1,1,1), 
		'encut': 400,
		'ismear': 0,
		'sigma': 0.1,
		'ediff': 1e-4,
		'prec': 'normal',
		'lcharg': True,
		'lwave': True,
		'nelmin': 4,
		'nelmdl': 6,
		'npar': 2,
		'algo': 'fast',
		'lreal': 'auto',
		'ispin': 2,
		'magmom': 0,
		'laechg': True,
		'istart': 0,
		},
	'In2O3_bader': {
		'xc': "PBE",
		'setups': {'In': '_d'},
		'kpts': (1,1,1),
		'encut': 400,
		'ismear': 0,
		'sigma': 0.05,
		'ediff': 1e-4,
		'prec': 'normal',
		'lcharg': True,
		'lwave': True,
		'nelmin': 4,
		'nelmdl': 6,
		'npar': 2,
		'algo': 'fast',
		'lreal': 'auto',
		'ispin': 2,
		'laechg': True,
		'istart': 0,
		},
	'ZnOCu_dimer': {
		'xc': "PBE",
		'kpts': (3,3,1),
		'encut': 400,
		'ismear': 0,
		'sigma': 0.1,
		'ediff': 1e-4,
		'lcharg': False,
		'lwave': True,
		'nelmin': 4,
		'nelmdl': 6,
		'npar': 2,
		'algo': 'fast',
		'lreal': 'auto',
		'ispin': 2,
		'nsw': 2000,
		'ediffg': -0.05,
		'iopt': 2,
		'ibrion': 3,
		'potim': 0,
		'ichain': 2,
		'drotmax': 6,
		'istart': 0,
		},
	'In2O3_dimer': {
		'xc': "PBE",
		'setups': {'In': '_d'},
		'kpts': (4,3,1),
		'encut': 400,
		'ismear': 0,
		'sigma': 0.05,
		'ediff': 1e-4,
		'prec': 'normal',
		'lcharg': False,
		'lwave': True,
		'nelmin': 4,
		'nelmdl': 6,
		'npar': 2,
		'algo': 'fast',
		'lreal': 'auto',
		'ispin': 2,
		'gamma': True,
		'nsw': 2000,
		'ediffg': -0.05,
		'iopt': 2,
		'ibrion': 3,
		'potim': 0,
		'ichain': 2,
		'drotmax': 6,
		'istart': 0,
		},
	'ZnOCu_NEB': {
		'images': 9,
		'xc': "PBE",
		'kpts': (3,3,1),
		'encut': 400,
		'ismear': 0,
		'sigma': 0.1,
		'ediff': 1e-4,
		'lcharg': False,
		'lwave': True,
		'nelmin': 4,
		'nelmdl': 6,
		'npar': 2,
		'algo': 'fast',
		'lreal': 'auto',
		'ispin': 2,
		'nsw': 2000,
		'ediffg': -0.10,
		'iopt': 1,
		'ibrion': 3,
		'potim': 0,
		'spring': -5,
		'lclimb': False,
		'istart': 0,
		'gamma': True,
		},
	'In2O3_NEB': {
		'images': 9,
		'setups': {'In': '_d'},
		'xc': "PBE",
		'kpts': (4,3,1),
		'encut': 400,
		'ismear': 0,
		'sigma': 0.05,
		'ediff': 1e-4,
		'lcharg': False,
		'lwave': True,
		'nelmin': 4,
		'nelmdl': 6,
		'npar': 2,
		'algo': 'fast',
		'lreal': 'auto',
		'ispin': 2,
		'nsw': 2000,
		'ediffg': -0.10,
		'iopt': 1,
		'ibrion': 3,
		'potim': 0,
		'spring': -5,
		'lclimb': False,
		'istart': 0,
		'gamma': True,
		},
	'In2O3_vib': {
		'xc': "PBE",
		'setups': {'In': '_d'},
		'kpts': (4,3,1),
		'encut': 400,
		'ismear': 0,
		'sigma': 0.05,
		'ediff': 1e-8,
		'prec': 'accurate',
		'lcharg': False,
		'lwave': False,
		'nelmin': 4,
		'nelmdl': 6,
		'npar': 2,
		'algo': 'fast',
		'lreal': 'auto',
		'ispin': 2,
		'gamma': True,
		'istart': 0,
		},
	'In2O3_dos': {
		'xc': "PBE",
		'setups': {'In': '_d'},
		'kpts': (4,3,1),
		'encut': 400,
		'ismear': 0,
		'sigma': 0.05,
		'ediff': 1e-4,
		'prec': 'normal',
		'lcharg': False,
		'lwave': False,
		'nelmin': 4,
		'nelmdl': 6,
		'npar': 2,
		'algo': 'fast',
		'lreal': 'auto',
		'ispin': 2,
		'gamma': True,
		'istart': 0,
		'lorbit': 12},
	'TiO2_step': {
		'xc': "PBE",
		'kpts': (2, 3, 1),
		'setups': {'Ti': '_sv'},
		#Parameters copied from RuO2 Step /home/work/ccei_biomass/users/jlym/SAC_project/RuO2/RuO2-Hs-step
		'lcharg': False,
		'lwave': False,
		'lvtot': False,
		'encut': 400,
		'algo': 'fast',
		'ismear': 0,
		'sigma': 0.1,
		'prec': 'normal',
		'lreal': 'auto',
		#'ropt': [2e-4, 2e-4, 2e-4],
		'istart': 0,
		'nelm': 400,
		'nelmdl': -10,
		'ediff': 1e-4,
		'ispin': 2,
		'isym': 0,
		'nsw': 1000,
		'isif': 2,
		'ibrion': 2,
		#'nfree': 2,
		'potim': 0.35,
		#ediffg uses VASP's internal optimizer. To update trajectory file this should be left out. It is controlled by LBFGS.run(fmax)
		#'ediffg': -0.05,
		'npar': 2,
		'lplane': True,
		#'nelmin': 4,
		'gamma': True,
		},
	'TiO2_anatase': {
		'xc': "PBE",
		'kpts': (2, 2, 1),
		'setups': {'Ti': '_sv'},
		#Parameters copied from RuO2 Step /home/work/ccei_biomass/users/jlym/SAC_project/RuO2/RuO2-Hs-step
		'lcharg': False,
		'lwave': False,
		'lvtot': False,
		'encut': 400,
		'algo': 'fast',
		'ismear': 0,
		'sigma': 0.1,
		'prec': 'normal',
		'lreal': 'auto',
		#'ropt': [2e-4, 2e-4, 2e-4],
		'istart': 0,
		'nelm': 400,
		'nelmdl': -10,
		'ediff': 1e-4,
		'ispin': 1,
		'isym': 0,
		'nsw': 1000,
		'isif': 2,
		'ibrion': 2,
		#'nfree': 2,
		'potim': 0.35,
		#ediffg uses VASP's internal optimizer. To update trajectory file this should be left out. It is controlled by LBFGS.run(fmax)
		#'ediffg': -0.05,
		'npar': 2,
		'lplane': True,
		#'nelmin': 4,
		'gamma': True,
		},
	'TiO2_anatase_spin': {
		'xc': "PBE",
		'kpts': (2, 2, 1),
		'setups': {'Ti': '_sv'},
		#Parameters copied from RuO2 Step /home/work/ccei_biomass/users/jlym/SAC_project/RuO2/RuO2-Hs-step
		'lcharg': False,
		'lwave': False,
		'lvtot': False,
		'encut': 400,
		'algo': 'fast',
		'ismear': 0,
		'sigma': 0.1,
		'prec': 'normal',
		'lreal': 'auto',
		#'ropt': [2e-4, 2e-4, 2e-4],
		'istart': 0,
		'nelm': 400,
		'nelmdl': -10,
		'ediff': 1e-4,
		'ispin': 2,
		'isym': 0,
		'nsw': 1000,
		'isif': 2,
		'ibrion': 2,
		#'nfree': 2,
		'potim': 0.35,
		#ediffg uses VASP's internal optimizer. To update trajectory file this should be left out. It is controlled by LBFGS.run(fmax)
		#'ediffg': -0.05,
		'npar': 2,
		'lplane': True,
		#'nelmin': 4,
		'gamma': True,
	},
	'TiO2_anatase_dipole': {
		'xc': "PBE",
		'kpts': (2, 2, 1),
		'setups': {'Ti': '_sv'},
		#Parameters copied from RuO2 Step /home/work/ccei_biomass/users/jlym/SAC_project/RuO2/RuO2-Hs-step
		'lcharg': False,
		'lwave': False,
		'lvtot': False,
		'encut': 400,
		'algo': 'very_fast',
		'ismear': 0,
		'sigma': 0.1,
		'prec': 'normal',
		'lreal': 'auto',
		'istart': 0,
		'nelm': 400,
		'nelmdl': -10,
		'ediff': 1e-4,
		'ispin': 1,
		'isym': 0,
		'nsw': 1000,
		'isif': 2,
		'ibrion': 2,
		#'nfree': 2,
		'potim': 0.35,
		#ediffg uses VASP's internal optimizer. To update trajectory file this should be left out. It is controlled by LBFGS.run(fmax)
		#'ediffg': -0.05,
		'npar': 2,
		'lplane': True,
		#'nelmin': 4,
		'gamma': True,
		#Dipole corrections
		'ldipol': True,
		'idipol': 3,
		'dipol': [0.5, 0.5, 0.5]
		},
	'FeO': {
		'xc': "PBE",
		'kpts': (4, 4, 1),
		'lcharg': False,
		'lwave': False,
		'lvtot': False,
		'encut': 400,
		'algo': 'fast',
		'ismear': 0,
		'sigma': 0.1,
		'prec': 'normal',
		'lreal': 'auto',
		'istart': 0,
		'nelm': 400,
		'nelmdl': -10,
		'ediff': 1e-4,
		'ispin': 2,
		'isym': 0,
		'nsw': 1000,
		'isif': 2,
		'ibrion': 2,
		'potim': 0.35,
		'npar': 2,
		'lplane': True,
		'gamma': True,
		'lorbit': 12,		
	},
	'neb': {
		'images': 9,
		'lwave': True,
		'iopt': 1,
		'ibrion': 3,
		'potim': 0,
		'spring': -5,
		'spring2': -5,
		'spower': 1,
		'ltangent': True,
		'lclimb': False,
		'ediffg': -0.15,
		'ichain': 0,
		'ispring': 1
		},
	'dimer':{'ediffg': -0.05,
		'iopt': 2,
		'ibrion': 3,
		'potim': 0,
		'ichain': 2,
		'drotmax': 6
		},
	'+U': {
		'ldau': True,
		'ldautype': 2,
        'lmaxmix': 4,
		},
	'bader': {
		'laechg': True,
		'lcharg': True,
		'lwave': True,
		},
	'Pt': {
		'xc': "PBE",
		'kpts': (3,3,1),
		'encut': 400,
		'ismear': 0,
		'sigma': 0.1,
		'ediff': 1e-4,
		'prec': 'normal',
		'gamma': True,
		'lcharg': False,
		'lwave': False,
		'nelmin': 4,
		'nelmdl': 6,
		'npar': 2,
		'algo': 'fast',
		'lreal': 'auto',
		'ispin': 2,
		'istart': 0,
		},
}

def add_neb(vasp_param, atoms_first = None, atoms_last = None, efirst = None, elast = None, **kwargs):
	for key, value in calc_dict['neb'].items():
		vasp_param[key] = value
	if efirst is not None:
		vasp_param['efirst'] = efirst 
	elif atoms_first is not None:
		vasp_param['efirst'] = atoms_first.get_potential_energy()

	if elast is not None:
		vasp_param['elast'] = elast
	elif atoms_last is not None:
		vasp_param['elast'] = atoms_last.get_potential_energy()

	for key, value in kwargs.items():
		vasp_param[key] = value
	return vasp_param

def add_dimer(vasp_param, **kwargs):
	for key, value in calc_dict['dimer'].items():
		vasp_param[key] = value
	for key, value in kwargs.items():
		vasp_param[key] = value
	return vasp_param

def add_plus_u(vasp_param, u_values, atoms = None, **kwargs):
	for key, value in calc_dict['+U'].items():
		vasp_param[key] = value
	#Resets default values with whatever's fed to it
	for key, value in kwargs.items():
		vasp_param[key] = value

	if atoms is not None and type(u_values) is dict:
		u_values_list = []
		elements = get_unique_list(atoms.get_chemical_symbols())
		for element in elements:
			try:
				u_values_list.append(u_values[element])
			except KeyError:
				#If the value was not specified, assumed to be zero
				u_values_list.append(0.0)
	else:
		u_values_list = list(u_values)
	vasp_param['ldauu'] = u_values_list
	return vasp_param

def add_bader(vasp_param, **kwargs):
	for key, value in calc_dict['bader'].items():
		vasp_param[key] = value
	for key, value in kwargs.items():
		vasp_param[key] = value
	return vasp_param

def print_vasp_param(calc):
	"""
	Prints the parameters of the vasp calculator that are not None will helpful text if available.
	"""
	print('-'*20)
	print('VASP calculator parameters:')
	max_space = 20
	type_params = [calc.input_params, calc.bool_params, calc.int_params, calc.dict_params, calc.exp_params, calc.float_params, calc.list_params, calc.special_params, calc.string_params]
	for type_param in type_params:
		for key, val in type_param.items():
			if val is not None:
				#help_text = help_dict[key]
				n_space = max_space - len(key) - len(str(val))
				#Use a single space if parameter and value are too long
				if n_space <= 0:
					n_space = 1
				#Checks for help entry in help_txt
				try:
					help_dict[key]
				except KeyError:
					help_txt = 'No help entry found for parameter %s' % key
				else:
					#Check if there is a nested dictionary
					if type(help_dict[key]) is dict:
						try:
							help_txt = help_dict[key][val]
						except KeyError:
							help_txt = 'No help entry found for parameter %s, key entry %r' % (key, val)
					else:
						help_txt = help_dict[key]
				#Account for the quotes when calculating the space to leave for the help text
				if type(val) is str:
					n_space -= 2					
				print('%s: %r%s%s' % (key, val, ' '*n_space, help_txt))
	print('-'*20)
   
def assign_magmom(atoms_obj, ispin = None):
	#Reference: http://kitchingroup.cheme.cmu.edu/dft-book/dft.html#orgheadline8
	magmom_dict = {'H': 1.,
				   'O': 2.,
				   'Fe': 2.22,
				   'Co': 1.72,
				   'Ni': 0.61}
	magmoms = [0]*len(atoms_obj)		  
	for i, atom in enumerate(atoms_obj):
		if magmom_dict.get(atom.symbol) is not None:
			magmoms[i] = magmom_dict.get(atom.symbol)
	return magmoms

def handle_restart(calc, atoms):
	"""Uses the istart and ispin parameters in calc to determine whether magmom should be set. """
	if calc.int_params['ispin'] == 2:
		if calc.int_params['istart'] == 0:
			magmom = [0] * len(atoms)
			calc.set(magmom = magmom)
#Help dictionary used for print_vasp_param
help_dict = { 
	#Boolean parameters
	'lcorr': {True: 'Harris Foulkes corrections applied.',
			  False: 'Harris Foulkes corrections NOT applied. [VASP Default]'},
	'lclimb': {True: 'Climbing NEB turned ON. [VASP Default]',
			   False: 'Climbing NEB turned OFF'},
	'lcharg': {True: 'LCHARG will be printed. [VASP Default]',
			   False: 'LCHARG will NOT be printed.'},
	'lnebcell': {True: 'Solid State NEB turned ON. Used with ISIF=3 and IOPT=3.',
				 False: 'Solid State NEB turned OFF. [VASP Default]'},
	'kgamma': {True: 'K point grid centered at gamma point. [VASP Default]',
			   False: 'K point grid not centered at gamma point.'},
	'lhfcalc': {True: 'Hartree-Fock type calculations will be performed.',
				False: 'Hartree-Fock type calculations will NOT be performed. [VASP Default]'},
	'lbeefens': {True: 'BEE energy contributions will be included in OUTCAR',
				 False: 'BEE energy corrections will NOT be included in OUTCAR'},
	'lvdw_ewald': {True: 'Lattice summation will be computed in E disp expression by means of Ewald''s summation.',
				   False: 'Lattice summation will NOT be computed in E disp expression by means of Ewald''s summation. [VASP Default]'},
	'lwave': {True: 'LWAVE will be printed. [VASP Default]',
			  False: 'LWAVE will NOT be printed.'},
	'lasync': {True: 'Communication will be overlapped with calculations.',
			   False: 'Communication will NOT be overlapped with calculations. [VASP Default]'},
	'lspectral': {True: 'The spectral method will be used. [VASP Default if NOMEGA > 2]',
				  False: 'The spectral method will NOT be used.'},
	'lvdw': {True: 'DFT-D2 method of Grimmie will be used.',
			 False: 'DFT-D2 method of Grimmie will NOT be used. [VASP Default]'},
	'ltangentold': {True: 'Old central difference tangent turned ON.',
					False: 'Old central difference tangent  turned OFF. [VASP Default]'},
	'ldau': {True: 'L(S)DA+U turned ON.',
			 False: 'L(S)DA+U turned OFF. [Vasp Default]'},
	'lsepb': {True: 'Charge density will be calculated for each band separately and written to PARCHG.nb.',
			  False: 'Charge density will be merged for all selected bands and written to PARCHG.ALLB or PARCHG. [VASP Default]'},
	'lcalcpol': {True: 'Evaluation of the Berry phase expressions for the macroscopic electronic polarization turned ON.',
				 False: 'Evaluation of the Berry phase expressions for the macroscopic electronic polarization turned OFF.'},
	'ldneb': {True: 'Modified double nudging turned ON.',
			  False: 'Modified double nudging turned OFF. [VASP Default]'},
	'lsepk': {True: 'Charge density of every k-point will be written to the files PARCHG.nk.',
			  False: 'Charge density of all k-points will be merged to a single file. [VASP Default]'},
	'lasph': {True: 'Non-spherical contributions from the gradient corrections inside the PAW spheres will be included.',
			  False: 'Non-spherical contributions from the gradient corrections inside the PAW spheres will NOT be included.'},
	'lwannier90': {True: 'Interface between VASP and WANNIER90 turned ON.',
				   False: 'Interface between VASP and WANNIER90 turned OFF. [VASP Default]'},
	'addgrid': {True: 'Additional support grid will be used to evaluate augmentation charges.',
				False: 'Additional support grid will NOT be used to evaluate augmentation charges.'},
	'ldiag': {True: 'Subspace rotation will be performed. [VASP Default]',
			  False: 'Subspace rotation will NOT be performed.'},			   
	'lelf': {True: 'Electron localization function will be written in ELFCAR.',
			 False: 'Electron localization function will NOT be written. [VASP Default]'},
	'loptics': {True: 'The frequency dependent dielectric matrix will be calculated after the electronic ground state has been determined.',
				False: 'The frequency dependent dielectric matrix will NOT be calculated. [VASP Default]'},
	'lcalceps': {True: 'The ion-clamped static dielectric tensor, the Born effective charge tensors, and the ion-clamped piezoelectric tensor of the system will be calculated from the self-consistent reponse to a finite electric field.',
				 False: 'The ion-clamped static dielectric tensor, the Born effective charge tensor, and the ion-clamped piezoelectric tensor of the system will NOT be calculated. [VASP Default]'},
	'luse_vdw': {True: 'VdW-DF functional will be used.',
				 False: 'VdW-DF function will NOT be used. [VASP Default]'},
	'lglobal': {True: 'The NEB will be optimized globally. [VASP Default]',
				False: 'The NEB images will be optimized image-by-image.'},
	'ldipol': {True: 'Dipole-dipole interaction correction will be applied to the potential.',
			   False: 'Dipole-dipole interaction correction will NOT be applied to the potential. [VASP Default]'},
	'lsorbit': {True: 'Spin-orbit coupling is turned ON.',
				False: 'Spin-orbit coupling is turned OFF. [VASP Default]'},
	'lautoscale': {True: 'The inverse curvature for VTST LBFGS will be automatically calculated.',
				   False: 'The inverse curvature for VTST LBFGS will NOT be automatically calculated.'},
	'lvhar': {True: 'Hartree potential will be written to LOCPOT.',
			  False: 'Hartree potential will NOT be written to LOCPOT.[VASP Default]'},
	'lbeefbas': {True: 'All BEEs will be printed in OUTCAR.',
				 False: 'All BEEs will not be printed in OUTCAR.'},
	'lvtot': {True: 'Total local potential will be written to the LOCPOT file.',
			  False: 'Total loc potential will NOT be written to the LOCPOT file. [VASP Default]'},
	'lscalu': {True: 'The parallel LU decomposition will be used in the orthonormalization of the wave functions.',
			   False: 'The parallel LU decomposition will NOT be used in the orthonormalization of the wave functions. [VASP Default]'},
	'llineopt': {True: 'A force-based line minimizer will be used for translation.',
				 False: 'A force-based line minimizer will NOT be used for translation. [VASP Default]'},
	'lpard': {True: 'Partial (band or k-point decomposed) charge densities will be evaluated.',
			  False: 'Partial (band or k-point decomposed) charge densities will NOT be evaluated. [VASP Default]'},
	'lthomas': {True: 'The decomposition of the exchange operator into a short range and a long range part will be based on Thomas-Fermi screening.',
				False: 'The decomposition of the exchange operator into a short range and a long range part will NOT be based on Thomas-Fermi screening. [VASP Default]'},
	'lrpa': {True: 'Local field effect will be induced on the Hartree level only.',
			 False: 'Charges of the Hartree and the exchange correlation potential are included. [VASP Default]'},
	'lsol': {True: 'Solvent effects will be added.',
			 False: 'Solvent effects will NOT be added. [VASP Default]'},
	'lplane': {True: 'Data distribution in real space will be done plane wise. [VASP Default]',
			   False: 'Data distribution in real space will NOT be done plane wise.'},
	'lepsilon': {True: 'The static dielectric matrix, ion-clamped piezoelectric tensor and the Born effective charges will be determined using density functional pertubation theory.',
				 False: 'The static dielectric matrix, ion-clamped piezoelectric tensor and the Born effective charges will NOT be determined using density functional pertubation theory. [VASP Default]'},
	'lscalapack': {True: 'scaLAPACK will be used.',
				   False: 'scaLAPACK will not be used. [VASP Default]'},
	'laechg': {True: 'Core charge will be written to AECCAR0 and valence charge will be written to AECCAR2.',
			   False: 'Core charge and valence charge will NOT be written. [VASP Default]'},
			  
	#Dictionary parameters
	'ldau_luj': 'Dictionary with L(S)DA+U parameters.',
			  
	#Exp parameters
	'ediff': 'Stopping-criterion for electronic SC-loop. Default = 1e-4',
	'symprec': 'Determines how accurate the positions in the POSCAR file must be. Default = 1e-5',
	'ediffg': 'Stopping-criterion for ionic SC-loop. Default = EDIFF * 10',
	'fdstep': 'Finite difference step for IOPT = 1 or 2. Default = 5e-3',
			  
	#Float parameters
	'tau': 'Surface tension parameter in Vaspsol',
	'emax': 'Energy range for DOSCAR file.',
	'falphadec': 'Factor to decrease alpha.',
	'encut': 'Planewave cutoff',
	'vdw_s8': 'Damping parameter for Grimme''s DFT-D3 dispersion correction.',
	'vdw_a1': 'Damping parameter for Grimme''s DFT-D3 dispersion correction.',
	'enaug': 'Kinetic energy cut-off for the augmentation changes.',
	'timestep': 'Dynamic timestep for IOPT = 3 and IOPT = 7.',
	'vdw_sr': 'Scaling parameter for Grimme''s DFT-D2 and DFT-D3 and Tkatchenko and Scheffler''s DFT-TS dispersion correction.',
	'maxmove': 'Max step for translation for IOPT > 0.',
	'ebreak': 'Absolute stopping criterion for optimization of eigenvalues.',
	'nelect': 'Total number of electrons.',
	'bmix': 'Cutoff wave vector for Kerker mixing scheme. Default = 1.0',
	'sdr': 'Finite difference for setting up Lanczos matrix and step size when translating.',
	'invcurve': 'Initial curvature for LBFGS (IOPT = 1)',
	'vdw_d': 'Global damping parameter for Grimme''s DFT-D2 and Tkatchenko and Scheffler''s DFT-TS dispersion corrections.',
	'ddr': '(DdR) Dimer separation.',
	'hfscreen': 'Attribute to change from PBE0 to HSE',
	'ftimeinc': 'Factor to increase dt.',
	'dfnmax': '(DFNMax) Rotational force below which dimer rotation stops.',
	'pstress': 'Added to the stress tensor.',
	'sdalpha': 'Ratio between force and step size for IOPT = 4',
	'pomass': 'Mass of ions in am',
	'clz': 'Electron count for core level shift.',
	'spring': 'Spring constant for NEB.',
	'zab_vdw': 'vdW-DF parameter.',
	'kspacing': 'Determines the number of k-points if the KPOINTS file is not present. KSPACING is the smallest allowed spacing between k-points.',
	'param2': 'Exchange parameter.',
	'amix_mag': 'Linear mixing parameter for the magnetization density. Default = 1.6',
	'aexx': 'Fraction of exact/DFT exchange.',
	'encutfock': 'FFT grid in the HF related routines.',
	'aggac': 'Fraction of gradient correction to correlation.',
	'ftimemax': 'Max time step.',
	'aldac': 'Fraction of LDA correlation energy.',
	'jacobian': 'Weight of lattice to atomic motion.',
	'vdw_scaling': 'Global scaling parameter for Grimme''s D-2 dispersion correction.',
	'vdw_cnradius': 'Cutoff radius for calculating coordination number in Grimme''s DFT-D3 dispersion correction.',
	'drotmax': '(DRotMax) Number of rotation steps per translation step.',
	'zval': 'Ionic valence.',
	'cshift': 'Sets the (small) complex shift η in Kramers-Kronig transformation.',
	'vdw_s6': 'Damping parameter for Grimme''s DFT-D2 and DFT-D3 and Tkatchenko and Scheffler''s DFT-TS dispersion corrections.',
	'deper': 'Relative stopping criterion for optimization of eigenvalues.',
	'bmix_mag': 'Cutoff wave vector for Kerker mixing scheme for the magnetization density. Default = 1.0',
	'aggax': 'Fraction of gradient correlation to exchange.',
	'stol': 'Convergence ratio for minimum eigenvalue.',
	'encutgw': 'Energy cutoff for response function.',
	'amin': 'Minimal mixing parameter in Kerker''s initial approximation to the charge dielectric function used in the Broyden/Pulay mixing scheme. Default = min(0.1, AMIX, AMIX_MAG)',
	'emin': 'Energy range for DOSCAR file.',
	'falpha': 'Parameter for velocity damping.',
	'weimin': 'Maximum weight for a band to be considered empty.',
	'ftimedec': 'Factor to decrease dt.',
	'vdw_radius': 'Cutoff radius for Grimme''s DFT-D2 and DFT-D3 and Tkatchenko and Scheffler''s DFT-TS dispersion corrections.',
	'dfnmin': '(DFNMin) Rotational force below which dimer is not rotated',
	'vdw_a2': 'Damping parameter for Grimme''s DFT-D3 dispersion correction.',
	'amix': 'Linear mixing parameter.',
	'param1': 'Exchange parameter.',
	'time': 'Time step for IALGO = 5X. Default = 0.4',
	'eb_k': 'Solvent permitivity in Vaspsol',
	'sigma': 'Width of smearing in eV. Default= 0.2',
	'efield': 'Size of the applied electric field in eV/A',
	'potim': 'Time-step for ion-motion in fs',

	#String parameters
	'algo': {'normal': 'Electronic minimization algorithm: Blocked Davidson iteration scheme (IALGO= 38) [VASP Default]',
			 'fast': 'Electronic minimization algorithm: The Blocked Davdison iteration scheme will be used for the initial phase followed by the RMM-DIIS',
			 'very_fast': 'Electronic minimization algorithm: RMM-DIIS (IALGO = 48).',
			 'conjugate': 'Electronic minimization algorithm: All band simultaneous update of orbitals',
			 'all': 'Electronic minimization algorithm: All band simultaneous update of orbitals',
			 'damped': 'Electronic minimization algorithm: damped velocity friction algorithm',
			 'subrot': 'Electronic minimization algorithm: Subspace rotation or diagonalization in the sub-space spanned by the calculated NBANDS orbitals.',
			 'exact': 'Electronic minimization algorithm: Exact diagonalization',
			 'diag': 'Electronic minimization algorithm: Exact diagonalization',
			 'eigenval': 'Electronic minimization algorithm: Recalculate one electron energies, density of state and perform selected postprocessing using the current orbitals.',
			 'nothing': 'Electronic minimization algorithm: Recalculate density of states (eigenvalues from WAVECAR) or perform other selected postprocessing using the current orbitals and one electron energies.',
			 'none': 'Electronic minimization algorithm: Recalculate density of states (eigenvalues from WAVECAR) or perform other selected postprocessing using the current orbitals and one electron energies.'},
	'gga': {'91': 'Type of GGA: Perdew-Wang 91',
			'PE': 'Type of GGA: Perdew-Burke-Erzenhof',
			'LM': 'Type of GGA: Langreth-Mehl-Hu',
			'PB': 'Type of GGA: Perdew-Becke'},
	'metagga': {'TPSS': 'See functional discussed in Sun et al. Phys. Rev. B 84, 035117 (2011)',
				'RTPSS': 'See functional discussed in Sun et al. Phys. Rev. B 84, 035117 (2011)',
				'M06L': 'See functional discussed in Y.Zhao and D.G. Truhlar, J. Chem. Phys. 125, 194101 (2006).',
				'MBJ': 'Modified Becke_johnson exchange potential.'},
	'prec': 'Default set of parameters used.',
	'system': 'Name of system.',
	'tebeg': 'Initial temperature for ab-initio molecular dynamics.',
	'teend': 'Final temperature for ab-initio molecular dynamics.',
	'precfock': 'Controls the FFT grid for the exact exchange (Hartree-Fock) routines.',
	
	#Special parameters
	'lreal': {True: 'Projection will be done in reciprocal space.',
			  False: 'Projection done in real space.',
			  'auto': 'Projection done in real space. Fully automatic optimization of projection operators. No user interference required.',
			  'on': 'Projection done in real space. Projection operators are re-optimized.'},
	
	#Integer parameters
	'ialgo': {-1: 'Integer selecting algorithm. Performance test.',
			  5: 'Integer selecting algorithm. Conjugate gradient algorithm. Steepest descent.',
			  6: 'Integer selecting algorithm. Conjugate gradient algorithm. Conjugated gradient.',
			  7: 'Integer selecting algorithm. Conjugate gradient algorithm. Preconditioned steepest descent.',
			  8: 'Integer selecting algorithm. Conjugate gradient algorithm. Preconditioned conjugated gradient. [VASP4.4 and older default]',
			  38: 'Integer selecting algorithm. Kosugi algorithm. [VASP4.5, 4.6 and 5.2 Default]',
			  44: 'Integer selecting algorithm. RMM-DIIS. Steepest descent eigenvalue minimization.',
			  46: 'Integer selecting algorithm. RMM-DIIS. Residuum-minimization + preconditioning',
			  48: 'Integer selecting algorithm. RMM-DIIS. Preconditioned residuum-minimization',
			  53: 'Integer selecting algorithm. Treat total free energy as variational quantity and minimize the functional completely selfconsistently. Damped MD with damping term automatically determined by the given time step.',
			  54: 'Integer selecting algorithm. Treat total free energy as variational quantity and minimize the functional completely selfconsistently. Damped MD (velocity quench or quickmin)',
			  58: 'Integer selecting algorithm. Treat total free energy as variational quantity and minimize the functional completely selfconsistently. Preconditioned conjugate gradient.',
			  2: 'Integer selecting algorithm. Orbitals and one-electron energies are kept fixed. One electron occupancies and electronic density of states (DOS) are, however, recalculated.',
			  3: 'Integer selecting algorithm. Orbitals (one-electron wavefunctions) are kept fixed. One-electron energies, one electron occupancies, band structure energies, and the electronic density of states (DOS) are, as well as, the total energy are recalculated for the present Hamiltonian.',
			  4: 'Integer selecting algorithm. Orbitals are updated by applying a sub-space rotation, i.e. the Hamiltonian is evaluated in the space spanned by the orbitals (read from WAVECAR), and one diagonalization in this space is performed. No optimization outside the subspace spanned by the orbitals is performed.',
			  15: 'Integer selecting algorithm. Conjugate gradient algorithm.',
			  16: 'Integer selecting algorithm. Conjugate gradient algorithm.',
			  17: 'Integer selecting algorithm. Conjugate gradient algorithm.',
			  18: 'Integer selecting algorithm. Conjugate gradient algorithm.',
			  28: 'Integer selecting algorithm. Conjugate gradient algorithm. Subspace-diagonalization before conjugate gradient algorithm.',
			  90: 'Integer selecting algorithm. Exact Diagonalization'},
	'ibrion': {-1: 'Determines how ions are updated and moved. No update; ions are not moved, but NSW outer loops are performed.',
			   0: 'Determines how ions are updated and moved. Standard ab-initio molecular dynamics.',
			   1: 'Determines how ions are updated and moved. A Quasi-Newton (variable metric) algorithm is used to relax the ions into their instantaneous groundstates.',
			   2: 'Determines how ions are updated and moved. A conjugate gradient algorithm is used to relax the ions into their instantaneous groundstate.',
			   3: 'Determines how ions are updated and moved.',
			   5: 'Determines how ions are updated and moved.',
			   6: 'Determines how ions are updated and moved.',
			   7: 'Determines how ions are updated and moved. It determines the Hessian matrix using density functional perturbation theory. Does NOT apply symmetry.',
			   8: 'Determines how ions are updated and moved. It determines the Hessian matrix using density functional perturbation theory. Apply symmetry.',
			   44: 'Determines how ions are updated and moved. Switches on the transition state optimization by means of the improved dimer method of Heyden et al.'},				  
	'icharg': {0: 'Determines how to construct the ''initial'' charge density. Calculate charge density from initial wave functions.',
			   1: 'Determines how to construct the ''initial'' charge density. Read from charge density from file CHGCAR.',
			   2: 'Determines how to construct the ''initial'' charge density. Take superposition of atomic charge densities.',
			   10: 'Determines how to construct the ''initial'' charge density. Non-selfconsistent calculation.',
			   11: 'Determines how to construct the ''initial'' charge density. Eigenvalues (for band structure plots), eigenfunctions or the DOS for a given charge density taken from CHGCAR. CHGCAR can be obtained from a selfconsistent calculation for a smaller k-points set.',
			   12: 'Determines how to construct the ''initial'' charge density. Non-self consistent calculations for a superposition of atomic charge densities.'},
	'idipol': {1: 'Dipole moment will be calculated only parallel to the direction of the 1st lattice vector.',
			   2: 'Dipole moment will be calculated only parallel to the direction of the 2nd lattice vector.',
			   3: 'Dipole moment will be calculated only parallel to the direction of the 3rd lattice vector.',
			   4: 'Full dipole moment in all directions will be calculated.'},
	'images': 'Number of images for NEB calculation',
	'iniwav': {0: 'Take ''jellium wave functions'', this means simply: fill wavefunction arrays with plane waves of lowest kinetic energy = lowest eigenvectors for a constant potential (''jellium'')',
			   1: 'Fill wavefunction arrays with random numbers. Use whenever possible. [VASP Default]'},	 # initial electr wf. : 0-lowe 1-rand
	'isif': {0: 'Force calculated. Stress tensor will NOT be calculated. Ions will be related. No change in cell shape. No change in cell volume.',
			 1: 'Force calculated. Only total pressure will be included in stress tensor. Ions will be related. No change in cell shape. No change in cell volume.',	   # calculate stress and what to relax
			 2: 'Force calculated. Stress tensor will be calculated. Ions will be related. No change in cell shape. No change in cell volume.',
			 3: 'Force calculated. Stress tensor will be calculated. Ions will be related. Allows change in cell shape. Allows change in cell volume.',
			 4: 'Force calculated. Stress tensor will be calculated. Ions will be related. Allows change in cell shape. No change in cell volume.',
			 5: 'Force calculated. Stress tensor will be calculated. Ions will NOT be related. Allows change in cell shape. No change in cell volume.',
			 6: 'Force calculated. Stress tensor will be calculated. Ions will NOT be related. Allows change in cell shape. Allows change in cell volume.',
			 7: 'Force calculated. Stress tensor will be calculated. Ions will NOT be related. No change in cell shape. Allows change in cell volume.'},
	'ismear': {-1: 'Fermi-smearing [VASP Default]',
			   0: 'Gaussian smearing',
			   1: 'Method of Methfessel-Paxton order 1',
			   -2: 'Partial occupancies are read in from INCAR, and kept fixed throughout run.',
			   -3: 'Make a loop over the different smearing-parameters supplied in the INCAR file.',
			   -4: 'Tetrahedron method without Blochl corrections.',
			   -5: 'Tetrahedron method without Blochl corrections.'},
	'ispin': {1: 'Non-spin polarized calculations. [VASP Default]',
			  2: 'Spin polarized calculations.'},
	'istart': {0: 'Start job ''from scratch''. Initialize wavefunction according to the flag INIWAV',
			   1: 'Continuation job. Read wavefunction from file WAVECAR.'},	 # startjob: 0-new 1-cont 2-samecut
	'isym': {-1: 'Symmetry turned OFF',
			 0: 'Symmetry turned OFF. Assumes that phi_k = phi*_-k and reduces the sampling of the Brillouin zone respectively.',
			 1: 'Symmetry turned ON.',
			 2: 'Symmetry turned ON. Memory conserving symmetrisation of the charge density is used.',
			 3: 'Symmetry turned ON. Charge density will be constructed by applying the relevant symmetry operations to the orbitals at the k-points in the irreducible part of the Brillouin zone.'},
#		'iwavpr',	 # prediction of wf.: 0-non 1-charg 2-wave 3-comb
#		'kpar',	   # k-point parallelization paramater
#		'ldauprint',  # 0-silent, 1-occ. matrix written to OUTCAR, 2-1+pot. matrix
#					  # written
	'ldautype': {1: 'L(S)DA+U applied. Liechtenstein approach',
				 2: 'L(S)DA+U applied. Dudarev approach',
				 4: 'L(S)DA+U applied. Liechtenstein (LDAU) approach'},   # L(S)DA+U: 1-Liechtenstein 2-Dudarev 4-Liechtenstein(LDAU)
#		'lmaxmix',	#
#		'lorbit',	 # create PROOUT
#		'maxmix',	 #
	'ngx': 'FFT mesh for wavefunctions, x',
	'ngxf': 'FFT mesh for charges x', 
	'ngy': 'FFT mesh for wavefunctions, y',
	'ngyf': 'FFT mesh for charges y',
	'ngz': 'FFT mesh for wavefunctions, z',
	'ngzf': 'FFT mesh for charges z',
	'nbands': 'Number of bands',
#		'nblk',	   # blocking for some BLAS calls (Sec. 6.5)
#		'nbmod',	  # specifies mode for partial charge calculation
	'nelm': 'Number of electronic steps (default 60)',
	'nelmdl': 'Number of initial electronic steps',
	'nelmin': 'Minimum number of initial electronic steps',
	'nfree': 'Number of steps per DOF when calculating Hessian using finite differences.',
#		'nkred',	  # define sub grid of q-points for HF with
#					  # nkredx=nkredy=nkredz
#		'nkredx',	  # define sub grid of q-points in x direction for HF
#		'nkredy',	  # define sub grid of q-points in y direction for HF
#		'nkredz',	  # define sub grid of q-points in z direction for HF
	'nomega': 'Number of frequency points',
	'nomegar': 'Number of frequency points on real axis',
	'npar': 'Parallelization over bands',
#		'nsim',	   # evaluate NSIM bands simultaneously if using RMM-DIIS
	'nsw': 'Number of steps for ionic updating.',
	'nupdown': 'Fix spin moment to specified value.',
#		'nwrite',	 # verbosity write-flag (how much is written)
#		'smass',	  # Nose mass-parameter (am)
#		'vdwgr',	  # extra keyword for Andris program
#		'vdwrn',	  # extra keyword for Andris program
#		'voskown',	# use Vosko, Wilk, Nusair interpolation
	# The next keywords pertain to the VTST add-ons from Graeme Henkelman's
	# group at UT Austin
#		'ichain',	 # Flag for controlling which method is being used (0=NEB,
#					  # 1=DynMat, 2=Dimer, 3=Lanczos) if ichain > 3, then both
#					  # IBRION and POTIM are automatically set in the INCAR file
#		'iopt',	   # Controls which optimizer to use.  for iopt > 0, ibrion = 3
#					  # and potim = 0.0
#		'snl',		# Maximum dimentionality of the Lanczos matrix
#		'lbfgsmem',   # Steps saved for inverse Hessian for IOPT = 1 (LBFGS)
#		'fnmin',	  # Max iter. before adjusting dt and alpha for IOPT = 7 (FIRE)
#		'icorelevel',  # core level shifts
#		'clnt',	   # species index
#		'cln',		# main quantum number of excited core electron
#		'cll',		# l quantum number of excited core electron
#		'ivdw',	   # Choose which dispersion correction method to use
	'nbandsgw': 'Number of bands for GW',
	'nbandso': 'Number of occupied bands for electron-hole treatment',
	'nbandsv': 'Number of virtual bands for electron-hole treatment',
	'ncore': 'Number of cores per band, equal to number of cores divided by npar',
#		'mdalgo',	 # Determines which MD method of Tomas Bucko to use
	'nedos': 'Number of grid points in DOS',
#		'turbo',	  # Ewald, 0 = Normal, 1 = PME

	#Input parameters
	'pp': 'Pseudopotential used.',
	'reciprocal': {True: 'Projection operators evaluated in reciprocal space.',
				   False: 'Projection operators evaluated in real space.'},
	#'kpts_nintersections'
	'setups': 'Different potentials used for different atoms.',
	'xc': 'Exchange functional',
	#'txt'
	'kpts': 'Number of k-points',
	'gamma': {True: 'K point will be placed at the gamma point.',
			  False: 'K point not necessarily placed at the gamma point.'}				  
}	
