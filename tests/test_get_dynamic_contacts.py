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
            self.assertEqual(len(lines), 24278)

            # Check that first two lines are comments
            self.assertEqual(lines[0][0], "#")
            self.assertEqual(lines[1][0], "#")
            self.assertEqual(lines[2][0], "0")

            # Check that 20 frames total are generated
            frames = set([int(l.split()[0]) for l in lines[2:]])
            self.assertEqual(frames, set(range(20)))

            # Check that there are no duplicate interactions, even if atom pairs are swapped
            contact_hashes = set()
            for line in lines[2:]:
                tokens = line.strip().split("\t")
                a1 = tokens[2]
                a2 = tokens[3]
                if a2 < a1:
                    a1, a2 = a2, a1
                chash = "\t".join(tokens[0:2] + [a1, a2])
                if chash in contact_hashes:
                    print(chash)
                self.assertFalse(chash in contact_hashes)
                contact_hashes.add(chash)

        os.remove(outfile)

    def test_empty_frames(self):
        """
        Tests that the multithreading doesnt discard a fragment when no contacts are found in the first frame
        (see issue #38).
        """
        outfile = "tests/5xnd_contacts.tsv"
        argv = ["--topology", "tests/5xnd_topology.pdb",
                "--trajectory", "tests/5xnd_trajectory.dcd",
                "--itypes", "sb",
                "--sele", "resname ASP LYS",
                "--output", outfile]
        get_dynamic_contacts.main(argv=argv)

        self.assertTrue(os.path.exists(outfile))
        with open(outfile) as f:
            lines = f.readlines()
            self.assertEqual(len(lines), 22)

            # Check that first two lines are comments
            self.assertEqual(lines[0][0], "#")
            self.assertEqual(lines[1][0], "#")
            self.assertEqual(lines[2][0], "1")

            # Check that 12 frames total are generated
            frames = set([int(l.split()[0]) for l in lines[2:]])
            self.assertEqual(frames, set([1, 2, 3, 4, 5, 8, 12, 13, 15, 16, 18, 19]))

        os.remove(outfile)


if __name__ == '__main__':
    unittest.main()
