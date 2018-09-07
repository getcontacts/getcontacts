#!/usr/bin/bash

# Example script for generating contact information. 

get_static_contacts.py \
  --structure 6cvo.pdb \
  --sele "nucleic" \
  --sele2 "nucleic" \
  --itypes hb \
  --output contacts.tsv \
  --ps_cutoff_dist 6.5 \
  --solv HOH
