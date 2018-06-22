__author__ = 'Rasmus Fonseca <fonseca.rasmus@gmail.com>'
__license__ = "Apache License 2.0"

import unittest
from contact_calc.contact_utils import *
import numpy as np
from numpy.linalg import norm


class TestContactUtils(unittest.TestCase):

    def test_points_to_vector(self):
        p1 = np.array([0, 0, 0])
        p2 = np.array([1, 0, 0])
        p3 = np.array([1, 1, 1])

        self.assertAlmostEqual(norm(points_to_vector(p1, p1)), 0)
        self.assertAlmostEqual(norm(points_to_vector(p1, p2)), 1)
        self.assertAlmostEqual(norm(points_to_vector(p2, p3)), math.sqrt(2))
        self.assertAlmostEqual(norm(points_to_vector(p3, p2) + points_to_vector(p2, p3)), 0)

    def test_vector_length(self):
        self.assertAlmostEqual(calc_vector_length(np.array([0, 0, 0])), 0)
        self.assertAlmostEqual(calc_vector_length(np.array([0, 3, 0])), 3)
        self.assertAlmostEqual(calc_vector_length(np.array([2, 0, 0])), 2)
        self.assertAlmostEqual(calc_vector_length(np.array([0, 0, 4])), 4)
        self.assertAlmostEqual(calc_vector_length(np.array([1, 1, 1])), math.sqrt(3))

    def test_calc_angle_between_vectors(self):
        self.assertAlmostEqual(calc_angle_between_vectors(np.array([1, 0, 0]), np.array([0, 1, 0])), 90, 5)
        self.assertAlmostEqual(calc_angle_between_vectors(np.array([1, 1, 0]), np.array([1, 1, 0])), 0, 5)
        self.assertAlmostEqual(calc_angle_between_vectors(np.array([1, 1, 0]), np.array([-1, -1, 0])), 180, 5)

    def test_calc_geom_distance(self):
        self.assertAlmostEqual(calc_geom_distance(np.array([1, 1, 1]), np.array([1, 1, 1])), 0, 5)
        self.assertAlmostEqual(calc_geom_distance(np.array([-1, -1, -1]), np.array([1, 1, 1])), math.sqrt(3)*2, 5)

if __name__ == '__main__':
    unittest.main()
