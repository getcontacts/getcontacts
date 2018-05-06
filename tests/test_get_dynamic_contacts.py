__author__ = 'Rasmus Fonseca <fonseca.rasmus@gmail.com>'
__license__ = "Apache License 2.0"

import unittest
import get_dynamic_contacts
import os


class TestGetDynamicContacts(unittest.TestCase):

    def test_5xnd(self):
        outfile = "tests/5xnd_contacts.tsv"
        argv = ("--topology tests/5xnd_topology.pdb "
                "--trajectory tests/5xnd_trajectory.dcd "
                "--output " + outfile + " "
                "--itypes all").split(" ")
        get_dynamic_contacts.main(argv=argv)

        self.assertTrue(os.path.exists(outfile))
        with open(outfile) as f:
            lines = f.readlines()
            self.assertEqual(len(lines), 22195)
            self.assertEqual(lines[0][0], "#")
            self.assertEqual(lines[1][0], "#")
            self.assertEqual(lines[2][0], "0")
            frames = set([int(l.split()[0]) for l in lines[2:]])
            self.assertEqual(frames, set(range(20)))

        os.remove(outfile)


if __name__ == '__main__':
    unittest.main()
