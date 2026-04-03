#!/usr/bin/bash

# Example script for generating contact information. 

get_static_contacts.py \
  --structure 6cvo.pdb \
  --sele "nucleic" \
  --sele2 "protein" \
  --ligand "resname AMP" \
  --itypes hb \
  --output contacts.tsv \
  --ps_cutoff_dist 6.5 \
  --hbond_cutoff_ang 70
