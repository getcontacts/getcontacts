__author__ = 'Rasmus Fonseca <fonseca.rasmus@gmail.com>'
__license__ = "Apache License 2.0"

import unittest
import get_dynamic_contacts
import get_contact_flare
import os


class TestGetDynamicContacts(unittest.TestCase):

    def test_5xnd_nolabel(self):
        contact_file = "tests/5xnd_contacts.tsv"
        argv = ("--topology tests/5xnd_topology.pdb "
                "--trajectory tests/5xnd_trajectory.dcd "
                "--output " + contact_file + " "
                "--itypes all").split(" ")
        get_dynamic_contacts.main(argv=argv)
        self.assertTrue(os.path.exists(contact_file))

        flare_file = "tests/5xnd_timeflare.tsv"
        argv = ("--input " + contact_file + " "
                "--output " + flare_file).split(" ")
        get_contact_flare.main(argv=argv)

        self.assertTrue(os.path.exists(flare_file))

        # TODO: test contents of flare_file

        os.remove(contact_file)
        os.remove(flare_file)


if __name__ == '__main__':
    unittest.main()
