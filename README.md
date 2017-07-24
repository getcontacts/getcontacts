# MDContactNetworks
Library for computing dynamic non-covalent contact networks in proteins throughout MD Simulation


## File format

MDContactNetworks is compatible with .nc and .dcd file formats for MD trajectories.


## Running the Code

The code is organized into two steps:

### 1. Computing non-covalent contacts in a protein throughout every frame of a MD Simulation fragment
   
   __Input:__ 

	Required Arguments:

	   TOP - Path to topology
	   TRAJ - Path to simulation trajectory fragment
	   OUTPUT_DIR - Path of output directory
	   INTERACTION_TYPES - User specifies what type of non-covalent interaction to compute using the following flags. 

		   -sb, Salt bridges
		   -pc, Pi-cation 
		   -ps, Pi-stacking
		   -ts, T-stacking
		   -vdw, Van Der Waals
		   -hb, Hydrogen Bonds
		   -lhb, Ligand Hydrogen Bonds

		   Hydrogen bonds are automatically stratified to following subtypes

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

		-cores <NUM_CORES> User can specify number of cpu cores to run computations on

		-solv <SOLVENT_ID> Specify solvent identifier in simulation other than the default of TIP3

		-chain <CHAIN_ID> Specify specific chain of protein to compute non-covalent contacts

		-ligand <LIGAND_NAME> Specify name of ligand if computing hydrogen bonds involving ligand

   
   __Output:__ List of non-covalent contacts of specified type for each frame of simulation fragment 

   __Examples:__

	python dynamic_contact_networks.py TOP.pdb TRAJ.nc -cores 12 -solv IP3 -chain A -ligand EJ4 -itype -hb -lhb


### 2. Visualizing dynamic contact networks

   __Input:__ Directory containing non-covalent contacts computed for each fragment of simulation

   
   __Output:__ JSON file that has keys representing interactions between a pair of atoms. Values being the time points that the interaction was formed in simulation. Formatted to be visualized with flareplots (https://gpcrviz.github.io/flareplot/)


