import unittest
import py_box3
import numpy as np

class TestConstants(unittest.TestCase):
	def test_base10_to_base3(self):
		np.testing.assert_array_equal(py_box3.base10_to_base3(n = 27, width = 4), np.array([1, 0, 0, 0]))
		np.testing.assert_array_equal(py_box3.base10_to_base3(n = 27, width = 6), np.array([0, 0, 1, 0, 0, 0]))

	def test_any_alpha(self):
		self.assertTrue(py_box3.any_alpha('123a'))
		self.assertTrue(py_box3.any_alpha('123A'))
		self.assertFalse(py_box3.any_alpha('123!@#_<>?/[]|\\'))

	def test_get_unique_list(self):
		self.assertListEqual(py_box3.get_unique_list(data = [1, 2, 2]), [1, 2])

	def test_get_RMSE(self):
		xs_data = [0.251193959,
				   0.690240061,
				   0.14934177,
				   0.711106899,
				   0.179392639,
				   0.481600941,
				   0.910675426,
				   0.089943402,
				   0.503349395,
				   0.530871672]

		xs_fit = [0.325905887,
				  0.716654719,
				  0.175646115,
				  0.7387606,
				  0.240873016,
				  0.555743172,
				  0.97102868,
				  0.146105241,
				  0.558958633,
				  0.559991138]
		
		self.assertAlmostEqual(py_box3.get_RMSE(xs_data = xs_data, xs_fit = xs_fit), 0.052678418)



if __name__ == '__main__':
	unittest.main()