
import unittest
import get_dynamic_contacts
import get_contact_ticc
import os


class TestGetContactTicc(unittest.TestCase):

    def test_5xnd_nolabel(self):
        """
        Take the 5xnd trajectory, compute contacts and then use TICC to extract two sets of frequencies
        """
        top_file = "tests/5xnd_topology.pdb"
        trj_file = "tests/5xnd_trajectory.dcd"
        contact_file = "tests/5xnd_contacts.tsv"
        args = "--topology %s --trajectory %s --output %s --itypes all" % (top_file, trj_file, contact_file)
        get_dynamic_contacts.main(argv=args.split(" "))
        self.assertTrue(os.path.exists(contact_file))

        ticc_output = (
            "tests/5xnd_segments.txt",
            "tests/5xnd_ticc"
        )
        args = "--input_contacts %s --clusters 2 --tab_output %s --frequency_output %s" % \
               ((contact_file,) + ticc_output)
        get_contact_ticc.main(argv=args.split(" "))
        ticc_output = (
            ticc_output[0],
            "tests/5xnd_ticc_resfreq_cluster000.tsv",
            "tests/5xnd_ticc_resfreq_cluster001.tsv"
        )

        for output_file in ticc_output:
            self.assertTrue(os.path.exists(output_file))

        for output_file in ticc_output:
            os.remove(output_file)
        os.remove(contact_file)


if __name__ == '__main__':
    unittest.main()


__author__ = 'Rasmus Fonseca <fonseca.rasmus@gmail.com>'
__license__ = "Apache License 2.0"
