__author__ = 'Rasmus Fonseca <fonseca.rasmus@gmail.com>'
__license__ = "Apache License 2.0"

import unittest
import get_static_contacts
import os


class TestGetStaticContacts(unittest.TestCase):

    def test_5y0t_interface(self):
        outfile = "tests/5Y0T_all.tsv"
        argv = ["--structure", "tests/5Y0T_h.pdb", "--output", outfile, "--itypes", "all",
                "--sele", "chain A", "--sele2", "chain B"]
        get_static_contacts.main(argv=argv)

        self.assertTrue(os.path.exists(outfile))
        with open(outfile) as f:
            lines = f.readlines()
            self.assertEqual(len(lines), 268)

            # Check that first two lines are comments
            self.assertEqual(lines[0][0], "#")
            self.assertEqual(lines[1][0], "#")
            self.assertEqual(lines[2][0], "0")

            # Check that the correct interaction types are present
            itypes = set([l.split()[1] for l in lines[2:]])
            self.assertEqual(itypes, set(["hbbb", "hbsb", "hbss", "hp", "lwb", "lwb2", "ps", "sb", "vdw", "wb", "wb2"]))

            # Check that all interactions span chain A and B
            chain_pairs = set([(l.split()[2][0], l.split()[3][0]) for l in lines[2:]])
            self.assertEqual(chain_pairs, set([("A", "B")]))

            # Check that there are no duplicate interactions, even if atom pairs are swapped
            contact_hashes = set()
            for line in lines[2:]:
                tokens = line.strip().split("\t")
                chash = tuple(sorted(tokens[1:]))
                if chash in contact_hashes:
                    print(chash)
                self.assertFalse(chash in contact_hashes)
                contact_hashes.add(chash)

        os.remove(outfile)


if __name__ == '__main__':
    unittest.main()
