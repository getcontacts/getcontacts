# MDContact
Library for computing dynamic non-covalent contact networks in proteins throughout MD Simulation.


## File format

MDContactNetworks is compatible with all topology and reimaged trajectory file formats readable in VMD.


## Running the Code

### 1. Computing non-covalent contacts in a protein throughout every frame of a MD Simulation fragment
   
   __Input:__ 

	Required Arguments:

	   	--topology TOPOLOGY - Path to topology
	   	--trajectory TRAJECTORY - Path to simulation trajectory fragment
	   	--output_dir OUTPUT_DIRECTORY - Path of output directory
	   	--itypes INTERACTION_TYPES - User specifies what type of non-covalent 
	   	interaction to compute using the following flags. 

		   -sb, Salt bridges
		   -pc, Pi-cation 
		   -ps, Pi-stacking
		   -ts, T-stacking
		   -vdw, Van Der Waals
		   -hb, Hydrogen Bonds
		   -lhb, Ligand Hydrogen Bonds

		   Hydrogen bonds are automatically stratified to following subtypes and 
		   written as output.

		   -hbbb, Backbone-backbone hydrogen bonds
		   -hbsb, Backbone-sidechain hydrogen bonds
		   -hbss, Sidechain-sidechain hydrogen bonds
		   -wb, Water-mediated hydrogen bond
		   -wb2, Extended water-mediated hydrogen bond
		   -hlb, Ligand-backbone hydrogen bonds
		   -hls, Ligand-sidechain hydrogen bonds
		   -lwb, Ligand water-mediated hydrogen bond
		   -lwb2, Ligand extended water-mediated hydrogen bond


	Optional Arguments:

		--cores NUM_CORES Number of cpu cores for parallelization,
		   default = 6

		--ligand LIGAND Resname of ligand molecule, default = None

		--sele SELECTION VMD selection query to compute contacts in specified 
		  region of protein, default = None

		--solv SOLVENT Solvent identifier in simulation, default = "TIP3"

	Arguments for adjusting geometric criteria:
		--sb_cutoff_dist SALT_BRIDGE_CUTOFF_DISTANCE
						cutoff for distance between anion and cation 
						atoms [default = 4.0 angstroms]
		--pc_cutoff_dist PI_CATION_CUTOFF_DISTANCE
						cutoff for distance between cation and centroid
						of aromatic ring [default = 6.0 angstroms]
		--pc_cutoff_ang PI_CATION_CUTOFF_ANGLE
						cutoff for angle between normal vector projecting
						from aromatic plane and vector from aromatic center
						to cation atom [default = 60 degrees]
		--ps_cutoff_dist PI_STACK_CUTOFF_DISTANCE
						cutoff for distance between centroids of two aromatic
						rings [default = 7.0 angstroms]
		--ps_cutoff_ang PI_STACK_CUTOFF_ANGLE
						cutoff for angle between the normal vectors projecting
						from each aromatic plane [default = 30 degrees]
		--ps_psi_ang PI_STACK_PSI_ANGLE
						cutoff for angle between normal vector projecting from
						aromatic plane 1 and vector between the two aromatic
						centroids [default = 45 degrees]
		--ts_cutoff_dist T_STACK_CUTOFF_DISTANCE
						cutoff for distance between centroids of two aromatic
						rings [default = 5.0 angstroms]
		--ts_cutoff_ang T_STACK_CUTOFF_ANGLE
						cutoff for angle between the normal vectors projecting
						from each aromatic plane minus 90 degrees 
						[default = 30 degrees]
		--ts_psi_ang T_STACK_PSI_ANGLE
						cutoff for angle between normal vector projecting from
						aromatic plane 1 and vector between the two aromatic
						centroids [default = 45 degrees]
		--hbond_cutoff_dist HBOND_CUTOFF_DISTANCE
						cutoff for distance between donor and acceptor atoms 
						[default = 3.5 angstroms]
		--hbond_cutoff_ang HBOND_CUTOFF_ANGLE
						cutoff for angle between donor hydrogen acceptor 
						[default = 70 degrees]
		--vdw_epsilon VDW_EPSILON
						amount of padding for calculating vanderwaals contacts 
						[default = 0.5 angstroms]

   
   __Output:__ Tables storing non-covalent contacts. Tab delimited rows are formatted to include 
   frame index, atom labels, and interaction types. 

   __Examples:__

	Salt bridges and hydrogen bonds for residues 100 to 160:
	python dynamic_contact_networks.py --topology TOP.pdb --trajectory TRAJ.nc --output_dir OUTPUT_DIR --cores 12 --solv IP3 --sele "chain A and resid 100 to 160" --ligand EJ4 --itype -sb -hb -lhb

	Pi-cation, pi-stacking, and vanderwaals contacts in the entire protein:
	python dynamic_contact_networks.py --topology TOP.psf --trajectory TRAJ.dcd --output_dir OUTPUT_DIR --cores 6 --itype -pc -ps -vdw

	Salt bridges and hydrogen bonds in the entire protein with modified distance cutoffs:
	python dynamic_contact_networks.py --topology TOP.mae --trajectory TRAJ.dcd --output_dir OUTPUT_DIR --cores 6 --sb_cutoff_dist 5.0 --hbond_cutoff_dist 4.5 --itype -sb -hb

