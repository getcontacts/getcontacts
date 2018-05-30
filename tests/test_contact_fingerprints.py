import unittest
import get_static_contacts
import get_contact_frequencies
import get_contact_fingerprints
import os


class TestGetContactFingerprints(unittest.TestCase):

    def test_5xnd_nolabel(self):
        """
        Take the four 5xnd models, compute contacts and then frequencies assuming that models 1 and 2 are from the same
        simulation and similarly for 3 and 4.
        """
        contact_io_files = [("tests/5xnd_%02d.pdb" % m, "tests/5xnd_%02d_contacts.tsv" % m) for m in range(1, 5)]
        for input_output in contact_io_files:
            args = "--structure %s --output %s --itypes all" % input_output
            get_static_contacts.main(argv=args.split(" "))
            self.assertTrue(os.path.exists(input_output[1]))

        resfreq_file_12 = "tests/5xnd_resfreq_12.tsv"
        args = "--input_files tests/5xnd_01_contacts.tsv tests/5xnd_02_contacts.tsv --output_file " + resfreq_file_12
        get_contact_frequencies.main(argv=args.split(" "))
        self.assertTrue(os.path.exists(resfreq_file_12))

        resfreq_file_34 = "tests/5xnd_resfreq_34.tsv"
        args = "--input_files tests/5xnd_03_contacts.tsv tests/5xnd_04_contacts.tsv --output_file " + resfreq_file_34
        get_contact_frequencies.main(argv=args.split(" "))
        self.assertTrue(os.path.exists(resfreq_file_12))

        fingerprint_files = (
            "tests/5xnd_fingerprint_table.tsv",
            "tests/5xnd_fingerprint_clustermap.png",
            "tests/5xnd_fingerprint_flare.json"
        )
        args = "--input_frequencies %s %s --table_output %s --plot_output %s --flare_output %s" % \
               ((resfreq_file_12, resfreq_file_34) + fingerprint_files)
        get_contact_fingerprints.main(argv=args.split(" "))
        for output_file in fingerprint_files:
            self.assertTrue(os.path.exists(output_file))

        for _, output_file in contact_io_files:
            os.remove(output_file)
        for output_file in fingerprint_files:
            os.remove(output_file)
        os.remove(resfreq_file_12)
        os.remove(resfreq_file_34)


if __name__ == '__main__':
    unittest.main()


__author__ = 'Rasmus Fonseca <fonseca.rasmus@gmail.com>'
__license__ = "Apache License 2.0"
