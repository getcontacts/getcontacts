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
	   OUTPUT - Path of output file 
	   INTERACTION_TYPE - User specifies what type of non-covalent interaction to compute using the following flags. 

		   -sb, Salt bridges
		   -pc, Pi-cation 
		   -ps, Pi-stacking
		   -ts, T-stacking
		   -vdw, Van Der Waals
		   -hbw, Hydrogen Bonds
		   -lhbw, Ligand Hydrogen Bonds

		   The following options require an additional -process flag

		   -hbbb, Backbone-backbone hydrogen bonds
		   -hbsb, Backbone-sidechain hydrogen bonds
		   -hbss, Sidechain-sidechain hydrogen bonds
		   -rw, Residue-water hydrogen bond
		   -wb, Water-mediated hydrogen bond
		   -wb2, Extended water-mediated hydrogen bond
		   -hlb, Ligand-backbone hydrogen bonds
		   -hls, Ligand-sidechain hydrogen bonds
		   -lw, Ligand-water hydrogen bond
		   -lwb, Ligand water-mediated hydrogen bond
		   -lwb2, Ligand extended water-mediated hydrogen bond


	Optional Arguments:

		-process <POST_PROCESS_FILE> To compute interaction types (hbbb, hbsb, etc) user must provide path to the already computed full hydrogen bonds interaction output hbw.txt

		-stride <STRIDE_VALUE> User can provide a stride value to subsample the trajectory

		-solv <SOLVENT_ID> Specify solvent identifier in simulation other than the default of TIP3

		-chain <CHAIN_ID> Specify specific chain of protein to compute non-covalent contacts

		-ligand <LIGAND_NAME> Specify name of ligand if computing hydrogen bonds involving ligand

   
   __Output:__ List of non-covalent contacts of specified type for each frame of simulation fragment 

   __Examples__


   1) python DynamicInteractionCalculator.py TOP.pdb TRAJ.nc hbw.txt -hbw -stride 5 -solv IP3 -chain A -ligand EJ4

   2) python DynamicInteractionCalculator.py TOP.pdb TRAJ.nc hbsb.txt -hbsb -process hbw.txt 


### 2. Visualizing dynamic contact networks

   __Input:__ Directory containing non-covalent contacts computed for each fragment of simulation

   
   __Output:__ JSON file that has keys representing interactions between a pair of atoms. Values being the time points that the interaction was formed in simulation. Formatted to be visualized with flareplots (https://gpcrviz.github.io/flareplot/)


