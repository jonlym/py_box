import unittest
from py_box3 import constants as c
import numpy as np

class TestConstants(unittest.TestCase):

	def test_R(self):
		self.assertEqual(c.R('J/mol/K'), 8.3144598)

	def test_h(self):
		self.assertEqual(c.h('J s', bar = False), 6.626070040e-34)
		self.assertEqual(c.h('J s', bar = True), 6.626070040e-34/(2.*np.pi))

	def test_convert_unit(self):
		self.assertEqual(c.convert_unit(num = 0., from_ = 'C', to = 'K'), 273.15)
		self.assertEqual(c.convert_unit(from_ = 'm', to = 'cm'), 100.)
		with self.assertRaises(ValueError):
			c.convert_unit(from_ = 'cm', to = 'J')

if __name__ == '__main__':
	unittest.main()