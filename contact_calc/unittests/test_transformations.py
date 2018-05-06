__author__ = 'Rasmus Fonseca <fonseca.rasmus@gmail.com>'
__license__ = "Apache License 2.0"

import unittest
import contact_calc.transformations as ct


class TestTransformations(unittest.TestCase):
    def setUp(self):
        self.input_lines = [
            "# total_frames:2 interaction_types:sb,vdw\n",
            "# Columns: frame, interaction_type, atom_1, atom_2[, atom_3[, atom_4]]\n",
            "0\tsb\tA:GLU:82:OE1\tA:ARG:76:NH2\n",
            "0\tsb\tA:GLU:82:OE2\tA:ARG:76:NH1\n",
            "1\tsb\tA:GLU:82:OE2\tA:ARG:76:NH2\n",
            "1\tvdw\tA:GLU:82:OE2\tA:ARG:76:NH2\n"
        ]

    def test_parse_contacts(self):
        contacts, num_frames = ct.parse_contacts(self.input_lines, set(['sb']))
        self.assertEqual(len(contacts), 3)
        self.assertEqual(num_frames, 2)
        self.assertTrue(all([c[1] == 'sb' for c in contacts]))
        self.assertEqual(type(contacts[0][0]), int)
        self.assertEqual(type(contacts[0][1]), str)
        self.assertEqual(type(contacts[0][2]), str)
        self.assertEqual(type(contacts[0][3]), str)

    def test_parse_residuelabels(self):
        pass

    def test_res_contacts(self):
        contacts, num_frames = ct.parse_contacts(self.input_lines, set(['sb']))
        rcontacts = ct.res_contacts(contacts)
        self.assertEqual(len(rcontacts), 2)
        self.assertEqual(rcontacts[0][0], 0)
        self.assertEqual(rcontacts[1][0], 1)
        self.assertEqual(rcontacts[1][1], "A:ARG:76")
        self.assertEqual(rcontacts[1][2], "A:GLU:82")


if __name__ == '__main__':
    unittest.main()
