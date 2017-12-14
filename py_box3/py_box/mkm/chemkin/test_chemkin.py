import unittest
import os
from py_box.mkm.chemkin import Chemkin

test_imp_path = os.path.join(os.environ['PYBOX'], 'py_box\\mkm\\chemkin\\test_INP')

class TestInputMethods(unittest.TestCase):

	def test_read_tube_inp(self):
		tube_dict = {'Ps': [1.0], 
					 'restart_max': 10, 
					 'sensitivity_analysis': False, 
					 'heat_transfer_coefficient': 0.0, 
					 'Ts': [973.0], 
					 'use_LSRs': False, 
					 'nnodes': 10, 
					 'use_omega': True, 
					 'use_coverage_effects': False, 
					 'isothermal': True, 
					 'reactor_type': 3, 
					 'save_transient': False, 
					 'reactant': "'NH3/GAS/'", 
					 'ttout': 0.01, 
					 'linear_T_ramp': 0, 
					 'SA_Vs': [650.0], 
					 'T_rise': 0.0, 
					 'set_equation_tolerance': False, 
					 'standard_T_and_P': True, 
					 'Qs': [1.77], 
					 'TPD_ramp': 2.0, 
					 'external_T': 923.0, 
					 'multi_input': False, 
					 'heat_transfer_area_to_volume': 3.571, 
					 'volume': 1.0, 
					 'use_BEPs': True, 
					 'use_iterative_solver': False, 
					 'ntdec': 10, 
					 'omega': 0.5, 
					 'upper_bandwidth': 0, 
					 'reaction_path_analysis_T': 900.0, 
					 'rtime': 10000000000.0, 
					 'non_negative_composition': True, 
					 'design_of_experiments': False, 
					 'use_binding_energy_corrections': False, 
					 'relative_tolerance': 1e-08, 
					 'MARI': "'N2/GAS/'", 
					 'reaction_path_analysis_mode': 1, 
					 'verbose_reaction_path_analysis': False, 
					 'use_different_activation_energy': False, 
					 'absolute_tolerance': 1e-10, 
					 'lower_bandwidth': 0, 
					 'n_runs': 1, 
					 'T_ref': 1.0}
		ck = Chemkin()
		self.assertDictEqual(ck.read_tube_inp(path = os.path.join(test_imp_path, 'tube.inp')), tube_dict)



if __name__ == '__main__':
	unittest.main()