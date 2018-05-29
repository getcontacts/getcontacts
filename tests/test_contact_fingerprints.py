__author__ = 'Rasmus Fonseca <fonseca.rasmus@gmail.com>'
__license__ = "Apache License 2.0"

import unittest
import get_static_contacts
import get_contact_frequencies
import os


class TestGetContactFingerprints(unittest.TestCase):

    def test_5xnd_nolabel(self):
        """
        Take the four 5xnd models, compute contacts and then frequencies assuming that models 1 and 2 are from the same
        simulation and similarly for 3 and 4.
        :return:
        """
        for model in range(1, 5):
            contact_file = "tests/5xnd_%02d_contacts.tsv" % model
            args = "--structure tests/5xnd_%02d.pdb --output %s --itypes all" % (model, contact_file)
            get_static_contacts.main(argv=args.split(" "))
            self.assertTrue(os.path.exists(contact_file))

        resfreq_file_12 = "tests/5xnd_resfreq_12.tsv"
        args = "--input_files tests/5xnd_01_contacts.tsv tests/5xnd_02_contacts.tsv --output_file " + resfreq_file_12
        get_contact_frequencies.main(argv=args.split(" "))
        self.assertTrue(os.path.exists(resfreq_file_12))

        resfreq_file_34 = "tests/5xnd_resfreq_34.tsv"
        args = "--input_files tests/5xnd_03_contacts.tsv tests/5xnd_04_contacts.tsv --output_file " + resfreq_file_34
        get_contact_frequencies.main(argv=args.split(" "))
        self.assertTrue(os.path.exists(resfreq_file_12))

        # TODO: Generate and test fingerprint

        # os.remove(contact_file)
        # os.remove(resfreq_file)


if __name__ == '__main__':
    unittest.main()
