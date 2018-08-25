#!/usr/bin/bash

# Example script for generating contact information. 

get_dynamic_contacts.py \
  --topology 5xnd_topology.pdb \
  --trajectory 5xnd_trajectory.dcd \
  --ligand "CA" \
  --itypes all \
  --output 5xnd_all-contacts.tsv

get_contact_trace.py \
  --input_contacts 5xnd_all-contacts.tsv \
  --interactions "A:CA:201.* A:GLU:60.*" \
                 "A:CA:201.* A:SER:56.*" \
                 "A:CA:201.* A:ASP:52.*" \
  --correlation_output 5xnd_CA_corr.png \
  --trace_output 5xnd_CA_trace.png
