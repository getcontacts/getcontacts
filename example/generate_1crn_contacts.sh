#!/usr/bin/bash

# Example script for generating contact information. 

pymol -c -d 'fetch 1crn, async=0; h_add; save 1crn_h.pdb'
rm -f 1crn.cif
python3 ../static_contacts.py --topology 1crn_h.pdb --all-interactions --output 1crn_all-contacts.tsv
