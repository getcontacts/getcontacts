
get_dynamic_contacts.py --topology dhfr_4m6l_top.psf --trajectory dhfr_4m6l_trj.dcd \
     --itypes hb --sele "protein or resname 21V NAP" --output dhfr_4m6l_hb_contacts.tsv

get_contact_frequencies.py --input_files dhfr_4m6l_hb_contacts.tsv \
     --itypes hbls hblb lwb lwb2 --label_file dhfr_4m6l_labels.tsv \
     --output dhfr_4m6l_freqs.tsv


get_dynamic_contacts.py --topology dhfr_8dfr_top.psf --trajectory dhfr_8dfr_trj.dcd \
     --itypes hb --sele "protein or resname NDP" --output dhfr_8dfr_hb_contacts.tsv

get_contact_frequencies.py --input_files dhfr_8dfr_hb_contacts.tsv \
     --itypes hbls hblb lwb lwb2 --label_file dhfr_8dfr_labels.tsv \
     --output dhfr_8dfr_freqs.tsv

get_contact_fingerprints.py --input_frequencies dhfr_{4m6l,8dfr}_freqs.tsv \
     --frequency_cutoff 0.6 --column_headers 4M6L 8DFR \
     --flare_output dhfr_compare.json \
     --plot_output dhfr_compare.png
