import unittest
import os
from py_box3.mkm.chemkin import Chemkin

test_imp_path = os.path.join(os.environ['PYBOX'], 'py_box\\mkm\\chemkin\\test_INP')

class TestInputMethods(unittest.TestCase):

	def test_read_tube_inp(self):
		ck = Chemkin()
		self.assertIsInstance(ck.read_tube_inp(path = os.path.join(test_imp_path, 'tube.inp')), dict)



if __name__ == '__main__':
	unittest.main()