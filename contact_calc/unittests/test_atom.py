__author__ = 'Rasmus Fonseca <fonseca.rasmus@gmail.com>'
__license__ = "Apache License 2.0"

import unittest
from contact_calc.atom import *


class TestContactUtils(unittest.TestCase):

    def test_init(self):
        a = Atom(1, "A", "ALA", 1, "N", "N")
        self.assertEqual(a.index, 1)
        self.assertEqual(a.chain, "A")
        self.assertEqual(a.resname, "ALA")
        self.assertEqual(a.resid, 1)
        self.assertEqual(a.name, "N")
        self.assertEqual(a.element, "N")
        self.assertAlmostEqual(a.vdwradius, 1.55)

        a = Atom(2, "A", "ALA", 1, "CA", "C")
        self.assertEqual(a.index, 2)
        self.assertEqual(a.name, "CA")
        self.assertEqual(a.element, "C")
        self.assertAlmostEqual(a.vdwradius, 1.7)

        a = Atom(1, "A", "ALA", 1, "CA", "X")
        self.assertEqual(a.name, "CA")
        self.assertEqual(a.element, "C")

        a = Atom(1, "A", "ALA", 1, "Gd", "X")
        self.assertEqual(a.name, "Gd")
        self.assertEqual(a.element, "GD")
        self.assertAlmostEqual(a.vdwradius, 2.79)

        a = Atom(1, "A", "ALA", 1, "XX", "X")
        self.assertEqual(a.name, "XX")
        self.assertEqual(a.element, "?")
        self.assertAlmostEqual(a.vdwradius, 1.7)


if __name__ == '__main__':
    unittest.main()
