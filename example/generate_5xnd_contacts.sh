#!/usr/bin/bash

# Example script for generating contact information. 

python3 ../contact_networks.py --topology 5xnd_topology.pdb --trajectory 5xnd_trajectory.dcd --all-interactions --output 5xnd_all-contacts.tsv
