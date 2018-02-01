import py_box3
import numpy as np
import unittest

class TestConstants(unittest.TestCase):

	def test_base10_to_base3(self):
		np.testing.assert_array_equal(py_box3.base10_to_base3(n = 10, width = 3), np.array([1., 0., 1.]))

	def test_any_alpha(self):
		self.assertFalse(py_box3.any_alpha('1234567890'))
		self.assertTrue(py_box3.any_alpha('1234567890F'))
		self.assertTrue(py_box3.any_alpha('1234567890f'))

	def test_get_RMSE(self):
		x = np.array(range(5))
		y = x + 1
		self.assertAlmostEqual(py_box3.get_RMSE(x, y), 1.)

	def test_spherical_to_xyz(self):
		r = 1.
		theta = 0.
		phi = 90.
		print(py_box3.spherical_to_xyz(r = r, theta = theta, phi = phi))
		np.testing.assert_allclose(py_box3.spherical_to_xyz(r = r, theta = theta, phi = phi), np.array([0.75, 0.433, 0.5]))	

if __name__ == '__main__':
	unittest.main()